! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
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
#define QUOTE(x) #x                                                                                 
#define PASTE(x,y) x ## y  

#ifndef INT
#define INT 4
#endif

#ifndef FLOAT
#define FLOAT 8
#endif

#include "prec.f90"
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
!> @details Usage:
!> @details   - choose material as hypela2
!> @details   - set statevariable 2 to index of homogenization
!> @details   - set statevariable 3 to index of microstructure
!> @details   - make sure the file "material.config" exists in the working directory
!> @details   - make sure the file "numerics.config" exists in the working directory
!> @details   - use nonsymmetric option for solver (e.g. direct profile or multifrontal sparse, 
!> @details     the latter seems to be faster!)
!> @details   - in case of ddm (domain decomposition) a SYMMETRIC solver has to be used, 
!> @details     i.e uncheck "non-symmetric"
!> @details  Marc subroutines used:
!> @details   - hypela2
!> @details   - plotv
!> @details   - quit
!> @details  Marc common blocks included:
!> @details   - concom: lovl, inc
!> @details   - creeps: timinc
!--------------------------------------------------------------------------------------------------
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
 character(len=*), parameter :: pathSep = achar(47)//achar(92)                                      ! forward and backward slash

 getSolverWorkingDirectoryName=''
 inputName=''
 inquire(5, name=inputName)                                                                         ! determine inputputfile
 getSolverWorkingDirectoryName=inputName(1:scan(inputName,pathSep,back=.true.))

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
 character(len=*), parameter :: pathSep = achar(47)//achar(92)                                      ! forward and backward slash
 integer(pInt) :: extPos

 getSolverJobName=''
 inputName=''
 inquire(5, name=inputName)                                                                         ! determine inputfile
 extPos = len_trim(inputName)-4
 getSolverJobName=inputName(scan(inputName,pathSep,back=.true.)+1:extPos)

end function getSolverJobName


end module DAMASK_interface

#include "IO.f90"
#include "libs.f90"
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
!> @details (1) F,R,U are only available for continuum and membrane elements (not for
!> @details     shells and beams).
!> @details
!> @details (2) Use the -> 'Plasticity,3' card(=update+finite+large disp+constant d)
!> @details     in the parameter section of input deck (updated Lagrangian formulation).
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
subroutine hypela2(d,g,e,de,s,t,dt,ngens,m,nn,kcus,matus,ndi,nshear,disp, &
                   dispt,coord,ffn,frotn,strechn,eigvn,ffn1,frotn1, &
                   strechn1,eigvn1,ncrd,itel,ndeg,ndm,nnode, &
                   jtype,lclass,ifr,ifu)
 use prec, only: &
   pReal, &
   pInt
 use numerics, only: &
!$ DAMASK_NumThreadsInt, &
   numerics_unitlength
 use FEsolving, only: &
   cycleCounter, &
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
 use math, only: &
   math_transpose33,&
   invnrmMandel
 use debug, only: &
   debug_level, &
   debug_LEVELBASIC, &
   debug_MARC, &
   debug_info, &
   debug_reset
 use mesh, only: &
   mesh_FEasCP, &
   mesh_element, &
   mesh_node0, &
   mesh_node, &
   mesh_Ncellnodes, &
   mesh_cellnode, &
   mesh_build_cellnodes, &
   mesh_build_ipCoordinates, &
   FE_Nnodes
 use CPFEM, only: &
   CPFEM_general, &
   CPFEM_init_done, &
   CPFEM_initAll, &
   CPFEM_CALCRESULTS, &
   CPFEM_AGERESULTS, &
   CPFEM_COLLECT, &
   CPFEM_RESTOREJACOBIAN, &
   CPFEM_BACKUPJACOBIAN

 implicit none
!$ include "omp_lib.h"                                                                              ! the openMP function library
 integer(pInt),                         intent(in) :: &                                             ! according to MSC.Marc 2012 Manual D
   ngens, &                                                                                         !< size of stress-strain law
   nn, &                                                                                            !< integration point number 
   ndi, &                                                                                           !< number of direct components
   nshear, &                                                                                        !< number of shear components
   ncrd, &                                                                                          !< number of coordinates
   itel, &                                                                                          !< dimension of F and R, either 2 or 3
   ndeg, &                                                                                          !< number of degrees of freedom
   ndm, &                                                                                           !< not specified in MSC.Marc 2012 Manual D
   nnode, &                                                                                         !< number of nodes per element
   jtype, &                                                                                         !< element type
   ifr, &                                                                                           !< set to 1 if R has been calculated
   ifu                                                                                              !< set to 1 if stretch has been calculated
 integer(pInt), dimension(2),           intent(in) :: &                                             ! according to MSC.Marc 2012 Manual D
   m, &                                                                                             !< (1) user element number, (2) internal element number
   matus, &                                                                                         !< (1) user material identification number, (2) internal material identification number
   kcus, &                                                                                          !< (1) layer number, (2) internal layer number 
   lclass                                                                                           !< (1) element class, (2) 0: displacement, 1: low order Herrmann, 2: high order Herrmann
 real(pReal),   dimension(*),           intent(in) :: &                                             ! has dimension(1) according to MSC.Marc 2012 Manual D, but according to example hypela2.f dimension(*)
   e, &                                                                                             !< total elastic strain
   de, &                                                                                            !< increment of strain
   dt                                                                                               !< increment of state variables
 real(pReal),   dimension(itel),        intent(in) :: &                                             ! according to MSC.Marc 2012 Manual D
   strechn, &                                                                                       !< square of principal stretch ratios, lambda(i) at t=n
   strechn1                                                                                         !< square of principal stretch ratios, lambda(i) at t=n+1
 real(pReal),   dimension(3,3),         intent(in) :: &                                             ! has dimension(itel,*) according to MSC.Marc 2012 Manual D, but we alway assume dimension(3,3)
   ffn, &                                                                                           !< deformation gradient at t=n
   ffn1                                                                                             !< deformation gradient at t=n+1
 real(pReal),   dimension(itel,*),      intent(in) :: &                                             ! according to MSC.Marc 2012 Manual D
   frotn, &                                                                                         !< rotation tensor at t=n
   eigvn, &                                                                                         !< i principal direction components for j eigenvalues at t=n
   frotn1, &                                                                                        !< rotation tensor at t=n+1
   eigvn1                                                                                           !< i principal direction components for j eigenvalues at t=n+1
 real(pReal),   dimension(ndeg,*),      intent(in) :: &                                             ! according to MSC.Marc 2012 Manual D
   disp, &                                                                                          !< incremental displacements
   dispt                                                                                            !< displacements at t=n (at assembly, lovl=4) and displacements at t=n+1 (at stress recovery, lovl=6)
 real(pReal),   dimension(ncrd,*),      intent(in) :: &                                             ! according to MSC.Marc 2012 Manual D
   coord                                                                                            !< coordinates
 real(pReal),   dimension(*),           intent(inout) :: &                                          ! according to MSC.Marc 2012 Manual D
   t                                                                                                !< state variables (comes in at t=n, must be updated to have state variables at t=n+1)
 real(pReal),   dimension(ndi+nshear),  intent(out) :: &                                            ! has dimension(*) according to MSC.Marc 2012 Manual D, but we need to loop over it
   s, &                                                                                             !< stress - should be updated by user
   g                                                                                                !< change in stress due to temperature effects
 real(pReal),   dimension(ngens,ngens), intent(out) :: &                                            ! according to MSC.Marc 2012 Manual D, but according to example hypela2.f dimension(ngens,*)
   d                                                                                                !< stress-strain law to be formed

!--------------------------------------------------------------------------------------------------
! Marc common blocks are in fixed format so they have to be reformated to free format (f90)
! Beware of changes in newer Marc versions

#include QUOTE(PASTE(../lib/MarcInclude/concom,Marc4DAMASK))                                        ! concom is needed for inc, lovl
#include QUOTE(PASTE(../lib/MarcInclude/creeps,Marc4DAMASK))                                        ! creeps is needed for timinc (time increment)

 logical :: cutBack
 real(pReal), dimension(6) ::   stress
 real(pReal), dimension(6,6) :: ddsdde
 integer(pInt) :: computationMode, i, cp_en, node, CPnodeID
 !$ integer(pInt) :: defaultNumThreadsInt                                                           !< default value set by Marc

 if (iand(debug_level(debug_MARC),debug_LEVELBASIC) /= 0_pInt) then
   write(6,'(a,/,i8,i8,i2)') ' MSC.MARC information on shape of element(2), IP:', m, nn
   write(6,'(a,2(1i))'), ' Jacobian:                      ', ngens,ngens 
   write(6,'(a,1i)'),    ' Direct stress:                 ', ndi  
   write(6,'(a,1i)'),    ' Shear stress:                  ', nshear  
   write(6,'(a,1i)'),    ' DoF:                           ', ndeg
   write(6,'(a,1i)'),    ' Coordinates:                   ', ncrd
   write(6,'(a,1i)'),    ' Nodes:                         ', nnode 
   write(6,'(a,1i)'),    ' Deformation gradient:          ', itel 
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' Deformation gradient at t=n:', &
                                 math_transpose33(ffn)
   write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' Deformation gradient at t=n+1:', &
                                 math_transpose33(ffn1)
 endif

 !$ defaultNumThreadsInt = omp_get_num_threads()                                                    ! remember number of threads set by Marc

 if (.not. CPFEM_init_done) call CPFEM_initAll(t(1),m(1),nn)

 !$ call omp_set_num_threads(DAMASK_NumThreadsInt)                                                  ! set number of threads for parallel execution set by DAMASK_NUM_THREADS

 computationMode = 0_pInt                                                                           ! save initialization value, since does not result in any calculation
 if (lovl == 4 ) then
   if(timinc < theDelta .and. theInc == inc ) computationMode = CPFEM_RESTOREJACOBIAN               ! first after cutback
 else                                                                                               ! stress requested (lovl == 6)
   cp_en = mesh_FEasCP('elem',m(1))
   if (cptim > theTime .or. inc /= theInc) then                                                     ! reached "convergence"
     terminallyIll = .false.
     cycleCounter = -1                                                                              ! first calc step increments this to cycle = 0
     if (inc == 0) then                                                                             ! >> start of analysis <<
       lastIncConverged = .false.                                                                   ! no Jacobian backup
       outdatedByNewInc = .false.                                                                   ! no aging of state
       lastMode = .false.                                                                           ! pretend last step was collection
       calcMode = .false.                                                                           ! pretend last step was collection
       !$OMP CRITICAL (write2out)
         write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> start of analysis..! ',m(1),nn
         flush(6)
       !$OMP END CRITICAL (write2out)
     else if (inc - theInc > 1) then                                                                ! >> restart of broken analysis <<
       lastIncConverged = .false.                                                                   ! no Jacobian backup
       outdatedByNewInc = .false.                                                                   ! no aging of state
       lastMode = .true.                                                                            ! pretend last step was calculation
       calcMode = .true.                                                                            ! pretend last step was calculation
       !$OMP CRITICAL (write2out)
         write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> restart of analysis..! ',m(1),nn
         flush(6)
       !$OMP END CRITICAL (write2out)
     else                                                                                           ! >> just the next inc <<
       lastIncConverged = .true.                                                                    ! request Jacobian backup
       outdatedByNewInc = .true.                                                                    ! request aging of state
       lastMode = .true.                                                                            ! assure last step was calculation
       calcMode = .true.                                                                            ! assure last step was calculation
       !$OMP CRITICAL (write2out)
         write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> new increment..! ',m(1),nn
         flush(6)
       !$OMP END CRITICAL (write2out)
     endif
   else if ( timinc < theDelta ) then                                                               ! >> cutBack <<
     terminallyIll = .false.
     cycleCounter = -1                                                                              ! first calc step increments this to cycle = 0
     calcMode = .true.                                                                              ! pretend last step was calculation
     !$OMP CRITICAL (write2out)
       write(6,'(a,i6,1x,i2)') '<< HYPELA2 >> cutback detected..! ',m(1),nn
       flush(6)
     !$OMP END CRITICAL (write2out)
   endif                                                                                            ! convergence treatment end

   calcMode(nn,cp_en) = .not. calcMode(nn,cp_en)                                                    ! ping pong (calc <--> collect)

   if ( calcMode(nn,cp_en) ) then                                                                   ! now --- CALC ---
     if ( lastMode /= calcMode(nn,cp_en) ) then                                                     ! first after ping pong
       call debug_reset()                                                                           ! resets debugging
       outdatedFFN1  = .false.
       cycleCounter  = cycleCounter + 1_pInt
       mesh_cellnode = mesh_build_cellnodes(mesh_node,mesh_Ncellnodes)                              ! update cell node coordinates
       call mesh_build_ipCoordinates()                                                              ! update ip coordinates
     endif
     if ( outdatedByNewInc ) then
       computationMode = ior(CPFEM_CALCRESULTS,CPFEM_AGERESULTS)
       outdatedByNewInc = .false.                                                                   ! reset flag
     else
       computationMode = CPFEM_CALCRESULTS
     endif
   else                                                                                             ! now --- COLLECT ---
     if ( lastMode /= calcMode(nn,cp_en) .and. & .not. terminallyIll) call debug_info()                                                                            ! first after ping pong reports (meaningful) debugging
     if ( lastIncConverged ) then
       computationMode = ior(CPFEM_COLLECT,CPFEM_BACKUPJACOBIAN)                                    ! collect and backup Jacobian after convergence
       lastIncConverged = .false.                                                                   ! reset flag
     else
       computationMode = CPFEM_COLLECT                                                              ! plain collect
     endif
     do node = 1,FE_Nnodes(mesh_element(2,cp_en))
       CPnodeID = mesh_element(4_pInt+node,cp_en)
       mesh_node(1:3,CPnodeID) = mesh_node0(1:3,CPnodeID) + numerics_unitlength * dispt(1:3,node)
     enddo
   endif

   theTime  = cptim                                                                                 ! record current starting time
   theDelta = timinc                                                                                ! record current time increment
   theInc   = inc                                                                                   ! record current increment number
   lastMode = calcMode(nn,cp_en)                                                                    ! record calculationMode
 endif

 call CPFEM_general(computationMode,ffn,ffn1,t(1),timinc,m(1),nn,stress,ddsdde)

!     Mandel: 11, 22, 33, SQRT(2)*12, SQRT(2)*23, SQRT(2)*13
!     Marc:   11, 22, 33, 12, 23, 13
!     Marc:   11, 22, 33, 12

 forall(i=1:ngens) d(1:ngens,i) = invnrmMandel(i)*ddsdde(1:ngens,i)*invnrmMandel(1:ngens)
 s(1:ndi+nshear) = stress(1:ndi+nshear)*invnrmMandel(1:ndi+nshear)
 g = 0.0_pReal
 if(symmetricSolver) d = 0.5_pReal*(d+transpose(d))
 
 !$ call omp_set_num_threads(defaultNumThreadsInt)                                                  ! reset number of threads to stored default value

end subroutine hypela2


!--------------------------------------------------------------------------------------------------
!> @brief sets user defined output variables for Marc
!> @details select a variable contour plotting (user subroutine).
!--------------------------------------------------------------------------------------------------
subroutine plotv(v,s,sp,etot,eplas,ecreep,t,m,nn,layer,ndi,nshear,jpltcd)
 use prec, only: &
   pReal, &
   pInt
 use mesh, only: &
   mesh_FEasCP
 use IO, only: &
   IO_error
 use homogenization, only: &
   materialpoint_results,&
   materialpoint_sizeResults
 
 implicit none
 integer(pInt),               intent(in) :: &
   m, &                                                                                             !< element number
   nn, &                                                                                            !< integration point number
   layer, &                                                                                         !< layer number
   ndi, &                                                                                           !< number of direct stress components
   nshear, &                                                                                        !< number of shear stress components
   jpltcd                                                                                           !< user variable index
 real(pReal),   dimension(*), intent(in) :: &
   s, &                                                                                             !< stress array
   sp, &                                                                                            !< stresses in preferred direction
   etot, &                                                                                          !< total strain (generalized) 
   eplas, &                                                                                         !< total plastic strain
   ecreep, &                                                                                        !< total creep strain
   t                                                                                                !< current temperature
 real(pReal),                 intent(out) :: &
   v                                                                                                !< variable

 if (jpltcd > materialpoint_sizeResults) call IO_error(700_pInt,jpltcd)                             ! complain about out of bounds error
 v = materialpoint_results(jpltcd,nn,mesh_FEasCP('elem', m))

end subroutine plotv
