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
!--------------------------------------------------------------------------------------------------
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief CPFEM engine
!--------------------------------------------------------------------------------------------------
module CPFEM
 use prec, only: &
   pReal, &
   pInt

 implicit none
 real(pReal), parameter ::                         CPFEM_odd_stress    = 1e15_pReal, &               !< return value for stress in case of ping pong dummy cycle
                                                   CPFEM_odd_jacobian  = 1e50_pReal                  !< return value for jacobian in case of ping pong dummy cycle

 real(pReal), dimension (:,:,:),   allocatable ::  CPFEM_cs                                          !< Cauchy stress
 real(pReal), dimension (:,:,:,:), allocatable ::  CPFEM_dcsdE                                       !< Cauchy stress tangent
 real(pReal), dimension (:,:,:,:), allocatable ::  CPFEM_dcsdE_knownGood                             !< known good tangent

 logical ::                                        CPFEM_init_done       = .false., &                !< remember whether init has been done already
                                                   CPFEM_init_inProgress = .false., &                !< remember whether first IP is currently performing init
                                                   CPFEM_calc_done       = .false.                   !< remember whether first IP has already calced the results
 logical, public, protected ::                     usePingPong           = .false.
 integer(pInt), parameter, public :: &
   CPFEM_CALCRESULTS     = 2_pInt**0_pInt, &
   CPFEM_AGERESULTS      = 2_pInt**1_pInt, &
   CPFEM_BACKUPJACOBIAN  = 2_pInt**2_pInt, &
   CPFEM_RESTOREJACOBIAN = 2_pInt**3_pInt, &
   CPFEM_COLLECT         = 2_pInt**4_pInt, &
   CPFEM_EXPLICIT        = 2_pInt**5_pInt

contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll(Temperature,element,IP)
 use prec, only: &
   prec_init
 use numerics, only: &
   numerics_init
 use debug, only: &
   debug_init
 use FEsolving, only: &
   FE_init
 use math, only: &
   math_init
 use mesh, only: &
   mesh_init
 use lattice, only: &
   lattice_init
 use material, only: &
   material_init
 use constitutive, only: &
   constitutive_init
 use crystallite, only: &
   crystallite_init
 use homogenization, only: &
   homogenization_init
 use IO, only: &
   IO_init
 use DAMASK_interface

 implicit none
 integer(pInt), intent(in) ::                        element, &                                    ! FE element number
                                                     IP                                            ! FE integration point number
 real(pReal), intent(in) ::                          Temperature                                   ! temperature
 real(pReal)   rnd
 integer(pInt) i,n

 ! initialization step (three dimensional stress state check missing?)

 if (.not. CPFEM_init_done) then
   call random_number(rnd)
   do i=1,int(256.0*rnd)
     n = n+1_pInt                                                                                  ! wasting random amount of time...
   enddo                                                                                           ! ...to break potential race in multithreading
   n = n+1_pInt
   if (.not. CPFEM_init_inProgress) then                                                           ! yes my thread won!
     CPFEM_init_inProgress = .true.
#ifdef Spectral
     call DAMASK_interface_init()                                                                  ! Spectral solver is interfacing to commandline
#endif
     call prec_init
     call IO_init
     call numerics_init
     call debug_init
     call math_init
     call FE_init
     call mesh_init(IP, element)                                                                   ! pass on coordinates to alter calcMode of first ip
     call lattice_init
     call material_init
     call constitutive_init
     call crystallite_init(Temperature)                                                            ! (have to) use temperature of first IP for whole model
     call homogenization_init(Temperature)
     call CPFEM_init
#ifndef Spectral
     call DAMASK_interface_init()                                                                  ! Spectral solver init is already done
#endif
     CPFEM_init_done = .true.
     CPFEM_init_inProgress = .false.
   else                                                                                            ! loser, loser...
     do while (CPFEM_init_inProgress)
     enddo
   endif
 endif

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init
 use, intrinsic :: iso_fortran_env                                                                 ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   pInt
 use IO, only: &
   IO_read_jobBinaryFile,&
   IO_read_jobBinaryIntFile, &
   IO_timeStamp, &
   IO_error
 use numerics, only: &
   DAMASK_NumThreadsInt
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive
 use FEsolving, only: &
   parallelExecution, &
   symmetricSolver, &
   restartRead, &
   modelName
 use mesh, only: &
   mesh_NcpElems, &
   mesh_Nelems, &
   mesh_maxNips
 use material, only: &
   homogenization_maxNgrains, &
   material_phase
 use constitutive, only: &
   constitutive_state0
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Lp0, &
   crystallite_dPdF0, &
   crystallite_Tstar0_v, &
   crystallite_localPlasticity
 use homogenization, only: &
   homogenization_sizeState, &
   homogenization_state0

 implicit none
 integer(pInt) i,j,k,l,m

 write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"

 if (any(.not. crystallite_localPlasticity) .and. (mesh_Nelems /= mesh_NcpElems)) call IO_error(600)
 if ((DAMASK_NumThreadsInt > 1_pInt)        .and. (mesh_Nelems /= mesh_NcpElems)) call IO_error(601)
 usePingPong = (any(.not. crystallite_localPlasticity) .or. (DAMASK_NumThreadsInt > 1_pInt))

 ! initialize stress and jacobian to zero 
 allocate(CPFEM_cs(6,mesh_maxNips,mesh_NcpElems)) ;                CPFEM_cs              = 0.0_pReal
 allocate(CPFEM_dcsdE(6,6,mesh_maxNips,mesh_NcpElems)) ;           CPFEM_dcsdE           = 0.0_pReal
 allocate(CPFEM_dcsdE_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dcsdE_knownGood = 0.0_pReal

 ! *** restore the last converged values of each essential variable from the binary file
 if (restartRead) then
   if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL (write2out)
      write(6,'(a)') '<< CPFEM >> restored state variables of last converged step from binary files'
     !$OMP END CRITICAL (write2out)
   endif

   call IO_read_jobBinaryIntFile(777,'recordedPhase',modelName,size(material_phase))
   read (777,rec=1) material_phase
   close (777)

   call IO_read_jobBinaryFile(777,'convergedF',modelName,size(crystallite_F0))
   read (777,rec=1) crystallite_F0
   close (777)

   call IO_read_jobBinaryFile(777,'convergedFp',modelName,size(crystallite_Fp0))
   read (777,rec=1) crystallite_Fp0
   close (777)

   call IO_read_jobBinaryFile(777,'convergedLp',modelName,size(crystallite_Lp0))
   read (777,rec=1) crystallite_Lp0
   close (777)

   call IO_read_jobBinaryFile(777,'convergeddPdF',modelName,size(crystallite_dPdF0))
   read (777,rec=1) crystallite_dPdF0
   close (777)

   call IO_read_jobBinaryFile(777,'convergedTstar',modelName,size(crystallite_Tstar0_v))
   read (777,rec=1) crystallite_Tstar0_v
   close (777)

   call IO_read_jobBinaryFile(777,'convergedStateConst',modelName)
   m = 0_pInt
   do i = 1,homogenization_maxNgrains; do j = 1,mesh_maxNips; do k = 1,mesh_NcpElems
     do l = 1,size(constitutive_state0(i,j,k)%p)
       m = m+1_pInt
       read(777,rec=m) constitutive_state0(i,j,k)%p(l)
     enddo
   enddo; enddo; enddo
   close (777)

   call IO_read_jobBinaryFile(777,'convergedStateHomog',modelName)
   m = 0_pInt
   do k = 1,mesh_NcpElems; do j = 1,mesh_maxNips
     do l = 1,homogenization_sizeState(j,k)
       m = m+1_pInt
       read(777,rec=m) homogenization_state0(j,k)%p(l)
     enddo
   enddo; enddo
   close (777)

   call IO_read_jobBinaryFile(777,'convergeddcsdE',modelName,size(CPFEM_dcsdE))
   read (777,rec=1) CPFEM_dcsdE
   close (777)
   restartRead = .false.
 endif

 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0) then
   write(6,'(a32,1x,6(i8,1x))') 'CPFEM_cs:              ', shape(CPFEM_cs)
   write(6,'(a32,1x,6(i8,1x))') 'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
   write(6,'(a32,1x,6(i8,1x))') 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
     write(6,*)
     write(6,*) 'parallelExecution:    ', parallelExecution
     write(6,*) 'symmetricSolver:      ', symmetricSolver
 endif
 flush(6)

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief perform initialization at first call, update variables and call the actual material model
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_general(mode, ffn, ffn1, Temperature, dt, element, IP, cauchyStress,&
     & jacobian, pstress, dPdF)
 ! note: cauchyStress = Cauchy stress cs(6) and jacobian = Consistent tangent dcs/dE
 use numerics, only:       defgradTolerance, &
                           iJacoStiffness
 use debug, only:          debug_level, &
                           debug_CPFEM, &
                           debug_levelBasic, &
                           debug_levelExtensive, &
                           debug_levelSelective, &
                           debug_e, &
                           debug_i, &
                           debug_stressMaxLocation, &
                           debug_stressMinLocation, &
                           debug_jacobianMaxLocation, &
                           debug_jacobianMinLocation, &
                           debug_stressMax, &
                           debug_stressMin, &
                           debug_jacobianMax, &
                           debug_jacobianMin
 use FEsolving, only:      parallelExecution, &
                           outdatedFFN1, &
                           terminallyIll, &
                           cycleCounter, &
                           theInc, &
                           theTime, &
                           theDelta, &
                           FEsolving_execElem, &
                           FEsolving_execIP, &
                           restartWrite
 use math, only:           math_identity2nd, &
                           math_mul33x33, &
                           math_det33, &
                           math_transpose33, &
                           math_I3, &
                           math_Mandel3333to66, &
                           math_Mandel66to3333, &
                           math_Mandel33to6, &
                           math_Mandel6to33
 use mesh, only:           mesh_FEasCP, &
                           mesh_NcpElems, &
                           mesh_maxNips, &
                           mesh_element, &
                           FE_Nips, &
                           FE_Nnodes, &
                           FE_geomtype
 use material, only:       homogenization_maxNgrains, &
                           microstructure_elemhomo, &
                           material_phase
 use constitutive, only:   constitutive_state0,constitutive_state
 use crystallite, only:    crystallite_partionedF,&
                           crystallite_F0, &
                           crystallite_Fp0, &
                           crystallite_Fp, &
                           crystallite_Lp0, &
                           crystallite_Lp, &
                           crystallite_dPdF0, &
                           crystallite_dPdF, &
                           crystallite_Tstar0_v, &
                           crystallite_Tstar_v
 use homogenization, only: homogenization_sizeState, &
                           homogenization_state, &
                           homogenization_state0, &
                           materialpoint_F, &
                           materialpoint_F0, &
                           materialpoint_P, &
                           materialpoint_dPdF, &
                           materialpoint_results, &
                           materialpoint_sizeResults, &
                           materialpoint_Temperature, &
                           materialpoint_stressAndItsTangent, &
                           materialpoint_postResults
 use IO, only:             IO_write_jobBinaryFile, &
                           IO_warning
 use DAMASK_interface
 
 implicit none
 integer(pInt), intent(in) ::                        element, &                                    ! FE element number
                                                     IP                                            ! FE integration point number
 real(pReal), intent(inout) ::                       Temperature                                   ! temperature
 real(pReal), intent(in) ::                          dt                                            ! time increment
 real(pReal), dimension (3,3), intent(in) ::         ffn, &                                        ! deformation gradient for t=t0
                                                     ffn1                                          ! deformation gradient for t=t1
 integer(pInt), intent(in) ::                        mode                                          ! computation mode  1: regular computation plus aging of results
 
 real(pReal), dimension(6), intent(out) ::           cauchyStress                                  ! stress vector in Mandel notation
 real(pReal), dimension(6,6), intent(out) ::         jacobian                                      ! jacobian in Mandel notation
 real(pReal), dimension (3,3), intent(out) ::        pstress                                       ! Piola-Kirchhoff stress in Matrix notation
 real(pReal), dimension (3,3,3,3), intent(out) ::    dPdF                                          ! 
 
 real(pReal)                                         J_inverse, &                                  ! inverse of Jacobian
                                                     rnd
 real(pReal), dimension (3,3) ::                     Kirchhoff, &                                  ! Piola-Kirchhoff stress in Matrix notation
                                                     cauchyStress33                                ! stress vector in Matrix notation
 real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                     H, &
                                                     jacobian3333                                  ! jacobian in Matrix notation
 integer(pInt)                                       cp_en, &                                      ! crystal plasticity element number
                                                     i, j, k, l, m, n                         
 logical                                             updateJaco                                    ! flag indicating if JAcobian has to be updated
 
 
 cp_en = mesh_FEasCP('elem',element)
 
 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt .and. cp_en == 1 .and. IP == 1) then
   !$OMP CRITICAL (write2out)
     write(6,'(/,a)') '#############################################'
     write(6,'(a1,a22,1x,f15.7,a6)') '#','theTime',        theTime,      '#'
     write(6,'(a1,a22,1x,f15.7,a6)') '#','theDelta',       theDelta,     '#'
     write(6,'(a1,a22,1x,i8,a13)')   '#','theInc',         theInc,       '#'
     write(6,'(a1,a22,1x,i8,a13)')   '#','cycleCounter',   cycleCounter, '#'
     write(6,'(a1,a22,1x,i8,a13)')   '#','computationMode',mode,         '#'
     write(6,'(a,/)') '#############################################'
     flush (6)
   !$OMP END CRITICAL (write2out)
 endif

     if (iand(mode, CPFEM_AGERESULTS) /= 0_pInt) then
       crystallite_F0  = crystallite_partionedF                                                    ! crystallite deformation (_subF is perturbed...)
       crystallite_Fp0 = crystallite_Fp                                                            ! crystallite plastic deformation
       crystallite_Lp0 = crystallite_Lp                                                            ! crystallite plastic velocity
       crystallite_dPdF0 = crystallite_dPdF                                                        ! crystallite stiffness
       crystallite_Tstar0_v = crystallite_Tstar_v                                                  ! crystallite 2nd Piola Kirchhoff stress 
       forall ( i = 1:homogenization_maxNgrains, &
                j = 1:mesh_maxNips, &
                k = 1:mesh_NcpElems ) &
         constitutive_state0(i,j,k)%p = constitutive_state(i,j,k)%p                                ! microstructure of crystallites
       if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
         !$OMP CRITICAL (write2out)
           write(6,'(a)') '<< CPFEM >> aging states'
           if (debug_e <= mesh_NcpElems .and. debug_i <= mesh_maxNips) then
             write(6,'(a,1x,i8,1x,i2,1x,i4,/,(12x,6(e20.8,1x)))') '<< CPFEM >> aged state of el ip grain',&
                                                             debug_e, debug_i, 1, constitutive_state(1,debug_i,debug_e)%p
             write(6,*)
           endif
         !$OMP END CRITICAL (write2out)
       endif
       !$OMP PARALLEL DO
         do k = 1,mesh_NcpElems
           do j = 1,mesh_maxNips
             if (homogenization_sizeState(j,k) > 0_pInt) &
               homogenization_state0(j,k)%p = homogenization_state(j,k)%p                          ! internal state of homogenization scheme
           enddo
         enddo
       !$OMP END PARALLEL DO

       ! * dump the last converged values of each essential variable to a binary file
       
       if (restartWrite) then 
         if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
          !$OMP CRITICAL (write2out)
            write(6,'(a)') '<< CPFEM >> writing state variables of last converged step to binary files'
          !$OMP END CRITICAL (write2out)
         endif
         
         call IO_write_jobBinaryFile(777,'recordedPhase',size(material_phase))
         write (777,rec=1) material_phase
         close (777)

         call IO_write_jobBinaryFile(777,'convergedF',size(crystallite_F0))
         write (777,rec=1) crystallite_F0
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergedFp',size(crystallite_Fp0))
         write (777,rec=1) crystallite_Fp0
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergedLp',size(crystallite_Lp0))
         write (777,rec=1) crystallite_Lp0
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergeddPdF',size(crystallite_dPdF0))
         write (777,rec=1) crystallite_dPdF0
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergedTstar',size(crystallite_Tstar0_v))
         write (777,rec=1) crystallite_Tstar0_v
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergedStateConst')
         m = 0_pInt
         do i = 1,homogenization_maxNgrains; do j = 1,mesh_maxNips; do k = 1,mesh_NcpElems
           do l = 1,size(constitutive_state0(i,j,k)%p)
             m = m+1_pInt
             write(777,rec=m) constitutive_state0(i,j,k)%p(l)
           enddo
         enddo; enddo; enddo
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergedStateHomog')
         m = 0_pInt
         do k = 1,mesh_NcpElems; do j = 1,mesh_maxNips
           do l = 1,homogenization_sizeState(j,k)
             m = m+1_pInt
             write(777,rec=m) homogenization_state0(j,k)%p(l)
           enddo
         enddo; enddo
         close (777)
         
         call IO_write_jobBinaryFile(777,'convergeddcsdE',size(CPFEM_dcsdE))
         write (777,rec=1) CPFEM_dcsdE
         close (777)
         
       endif
     endif

     if (iand(mode, CPFEM_EXPLICIT) /= 0_pInt) then                                                            ! Abaqus explicit skips collect
       materialpoint_Temperature(IP,cp_en) = Temperature
       materialpoint_F0(1:3,1:3,IP,cp_en) = ffn
       materialpoint_F(1:3,1:3,IP,cp_en) = ffn1
     endif

     if (iand(mode, CPFEM_CALCRESULTS) /= 0_pInt) then
     !*** deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
     
     if (terminallyIll .or. outdatedFFN1 .or. any(abs(ffn1 - materialpoint_F(1:3,1:3,IP,cp_en)) > defgradTolerance)) then
!        if (.not. terminallyIll .and. .not. outdatedFFN1) then 
       if (any(abs(ffn1 - materialpoint_F(1:3,1:3,IP,cp_en)) > defgradTolerance)) then 
         if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /=  0_pInt) then
           !$OMP CRITICAL (write2out)
             write(6,'(a,1x,i8,1x,i2)') '<< CPFEM >> OUTDATED at el ip',cp_en,IP
             write(6,'(a,/,3(12x,3(f10.6,1x),/))') '<< CPFEM >> FFN1 old:',math_transpose33(materialpoint_F(1:3,1:3,IP,cp_en))
             write(6,'(a,/,3(12x,3(f10.6,1x),/))') '<< CPFEM >> FFN1 now:',math_transpose33(ffn1)
           !$OMP END CRITICAL (write2out)
         endif
         outdatedFFN1 = .true.
       endif
       call random_number(rnd)
       if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
       CPFEM_cs(1:6,IP,cp_en) = rnd*CPFEM_odd_stress
       CPFEM_dcsde(1:6,1:6,IP,cp_en) = CPFEM_odd_jacobian*math_identity2nd(6)
     
     
     !*** deformation gradient is not outdated
     
     else
       updateJaco = mod(cycleCounter,iJacoStiffness) == 0
       
       !* no parallel computation, so we use just one single element and IP for computation
       
       if (.not. parallelExecution) then
         FEsolving_execElem(1)     = cp_en
         FEsolving_execElem(2)     = cp_en
         if (.not. microstructure_elemhomo(mesh_element(4,cp_en)) .or. &                           ! calculate unless homogeneous
                  (microstructure_elemhomo(mesh_element(4,cp_en)) .and. IP == 1_pInt)) then        ! and then only first IP
           FEsolving_execIP(1,cp_en) = IP
           FEsolving_execIP(2,cp_en) = IP
           if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /=  0_pInt) then
             !$OMP CRITICAL (write2out)
               write(6,'(a,i8,1x,i2)') '<< CPFEM >> calculation for el ip ',cp_en,IP
             !$OMP END CRITICAL (write2out)
           endif
           call materialpoint_stressAndItsTangent(updateJaco, dt)                                    ! calculate stress and its tangent
           call materialpoint_postResults(dt)                                                        ! post results
         endif
         
       !* parallel computation and calulation not yet done
       
       elseif (.not. CPFEM_calc_done) then
         if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
           !$OMP CRITICAL (write2out)
             write(6,'(a,i8,a,i8)') '<< CPFEM >> calculation for elements ',FEsolving_execElem(1),' to ',FEsolving_execElem(2)
           !$OMP END CRITICAL (write2out)
         endif
         call materialpoint_stressAndItsTangent(updateJaco, dt)                                    ! calculate stress and its tangent (parallel execution inside)
         call materialpoint_postResults(dt)                                                        ! post results
         CPFEM_calc_done = .true.
       endif
       
               !*** map stress and stiffness (or return odd values if terminally ill)
       
       if ( terminallyIll ) then
         call random_number(rnd)
         if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
         CPFEM_cs(1:6,IP,cp_en) = rnd * CPFEM_odd_stress
         CPFEM_dcsde(1:6,1:6,IP,cp_en) = CPFEM_odd_jacobian * math_identity2nd(6)
       else  
         if (microstructure_elemhomo(mesh_element(4,cp_en)) .and. IP > 1_pInt) then                ! me homogenous? --> copy from first IP
           materialpoint_P(1:3,1:3,IP,cp_en) = materialpoint_P(1:3,1:3,1,cp_en)
           materialpoint_F(1:3,1:3,IP,cp_en) = materialpoint_F(1:3,1:3,1,cp_en)
           materialpoint_dPdF(1:3,1:3,1:3,1:3,IP,cp_en) = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,cp_en)
           materialpoint_results(1:materialpoint_sizeResults,IP,cp_en) = materialpoint_results(1:materialpoint_sizeResults,1,cp_en)
         endif
         ! translate from P to CS
         Kirchhoff = math_mul33x33(materialpoint_P(1:3,1:3,IP,cp_en), math_transpose33(materialpoint_F(1:3,1:3,IP,cp_en)))
         J_inverse  = 1.0_pReal / math_det33(materialpoint_F(1:3,1:3,IP,cp_en))
         CPFEM_cs(1:6,IP,cp_en) = math_Mandel33to6(J_inverse * Kirchhoff)

       !  translate from dP/dF to dCS/dE
         H = 0.0_pReal
         do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
           H(i,j,k,l) = H(i,j,k,l) + &
                         materialpoint_F(j,m,IP,cp_en) * &
                         materialpoint_F(l,n,IP,cp_en) * &
                         materialpoint_dPdF(i,m,k,n,IP,cp_en) - &
                         math_I3(j,l) * materialpoint_F(i,m,IP,cp_en) * materialpoint_P(k,m,IP,cp_en) + &
                         0.5_pReal * (math_I3(i,k) * Kirchhoff(j,l) + math_I3(j,l) * Kirchhoff(i,k) + &
                                    math_I3(i,l) * Kirchhoff(j,k) + math_I3(j,k) * Kirchhoff(i,l))
         enddo; enddo; enddo; enddo; enddo; enddo
         do i=1,3; do j=1,3; do k=1,3; do l=1,3
           H_sym(i,j,k,l) = 0.25_pReal * (H(i,j,k,l) + H(j,i,k,l) + H(i,j,l,k) + H(j,i,l,k))
         enddo; enddo; enddo; enddo
         CPFEM_dcsde(1:6,1:6,IP,cp_en) = math_Mandel3333to66(J_inverse * H_sym)
       endif
     endif
   endif

   ! --+>> COLLECTION OF FEM INPUT WITH RETURNING OF RANDOMIZED ODD STRESS AND JACOBIAN <<+-- 


     if (iand(mode, CPFEM_BACKUPJACOBIAN) /= 0_pInt) &
       CPFEM_dcsde_knownGood = CPFEM_dcsde
     if (iand(mode, CPFEM_RESTOREJACOBIAN) /= 0_pInt) &
       CPFEM_dcsde = CPFEM_dcsde_knownGood
     
     if (iand(mode, CPFEM_COLLECT) /= 0_pInt) then
       call random_number(rnd)
       if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
       materialpoint_Temperature(IP,cp_en) = Temperature
       materialpoint_F0(1:3,1:3,IP,cp_en) = ffn
       materialpoint_F(1:3,1:3,IP,cp_en) = ffn1
       CPFEM_cs(1:6,IP,cp_en) = rnd * CPFEM_odd_stress
       CPFEM_dcsde(1:6,1:6,IP,cp_en) = CPFEM_odd_jacobian * math_identity2nd(6)
       CPFEM_calc_done = .false.
     endif

 
 cauchyStress = CPFEM_cs(1:6,IP,cp_en)
 jacobian = CPFEM_dcsdE(1:6,1:6,IP,cp_en)
 pstress = materialpoint_P(1:3,1:3,IP,cp_en)
 dPdF = materialpoint_dPdF(1:3,1:3,1:3,1:3,IP,cp_en)
 if (theTime > 0.0_pReal) then
   Temperature = materialpoint_Temperature(IP,cp_en)  ! homogenized result except for potentially non-isothermal starting condition.
 endif

 if (mode < 3 .and. iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt &
                    .and. ((debug_e == cp_en .and. debug_i == IP) &
                            .or. .not. iand(debug_level(debug_CPFEM), debug_levelSelective) /= 0_pInt)) then
   !$OMP CRITICAL (write2out)
     write(6,'(a,i8,1x,i2,/,12x,6(f10.3,1x)/)') '<< CPFEM >> stress/MPa at el ip ', cp_en, IP, cauchyStress/1.0e6_pReal
     write(6,'(a,i8,1x,i2,/,6(12x,6(f10.3,1x)/))') '<< CPFEM >> Jacobian/GPa at el ip ', cp_en, IP, transpose(jacobian)/1.0e9_pReal
     flush(6)
   !$OMP END CRITICAL (write2out)
 endif
 
 
 !*** warn if stiffness close to zero
 
 if (all(abs(jacobian) < 1e-10_pReal)) then
   call IO_warning(601,cp_en,IP)
 endif
 
 
 !*** remember extreme values of stress and jacobian
 
 if (mode < 3) then
   cauchyStress33 = math_Mandel6to33(cauchyStress)
   if (maxval(cauchyStress33) > debug_stressMax) then                        
     debug_stressMaxLocation = [cp_en, IP]
     debug_stressMax = maxval(cauchyStress33)
   endif
   if (minval(cauchyStress33) < debug_stressMin) then
     debug_stressMinLocation = [cp_en, IP]
     debug_stressMin = minval(cauchyStress33)
   endif
   jacobian3333 = math_Mandel66to3333(jacobian)
   if (maxval(jacobian3333) > debug_jacobianMax) then
     debug_jacobianMaxLocation = [cp_en, IP]
     debug_jacobianMax = maxval(jacobian3333)
   endif
   if (minval(jacobian3333) < debug_jacobianMin) then
     debug_jacobianMinLocation = [cp_en, IP]
     debug_jacobianMin = minval(jacobian3333)
   endif
 endif
 
end subroutine CPFEM_general

end module CPFEM
