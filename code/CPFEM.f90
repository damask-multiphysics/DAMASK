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
 private
#if defined(Marc4DAMASK) || defined(Abaqus)
 real(pReal),                      parameter,   private :: &
   CPFEM_odd_stress    = 1e15_pReal, &                                                              !< return value for stress in case of ping pong dummy cycle
   CPFEM_odd_jacobian  = 1e50_pReal                                                                 !< return value for jacobian in case of ping pong dummy cycle
 real(pReal), dimension (:,:,:),   allocatable, private :: &
   CPFEM_cs                                                                                         !< Cauchy stress
 real(pReal), dimension (:,:,:,:), allocatable, private :: &
   CPFEM_dcsdE                                                                                      !< Cauchy stress tangent
 real(pReal), dimension (:,:,:,:), allocatable, private :: &
   CPFEM_dcsdE_knownGood                                                                            !< known good tangent
#endif
 integer(pInt),                                 public :: &
   cycleCounter =  0_pInt, &                                                                        !< needs description
   theInc       = -1_pInt, &                                                                        !< needs description
   lastLovl     =  0_pInt, &                                                                        !< lovl in previous call to marc hypela2
   lastStep     =  0_pInt                                                                           !< kstep in previous call to abaqus umat
 real(pReal),                                   public :: &
   theTime      = 0.0_pReal, &                                                                      !< needs description
   theDelta     = 0.0_pReal    
 logical,                                       public :: & 
   outdatedFFN1      = .false., &                                                                   !< needs description
   lastIncConverged  = .false., &                                                                   !< needs description
   outdatedByNewInc  = .false.                                                                      !< needs description

 logical,                                       public, protected :: &
   CPFEM_init_done       = .false.                                                                  !< remember whether init has been done already
 logical,                                       private :: &
   CPFEM_calc_done       = .false.                                                                  !< remember whether first ip has already calced the results

 integer(pInt), parameter,                      public :: &
   CPFEM_COLLECT         = 2_pInt**0_pInt, &
   CPFEM_CALCRESULTS     = 2_pInt**1_pInt, &
   CPFEM_AGERESULTS      = 2_pInt**2_pInt, &
   CPFEM_BACKUPJACOBIAN  = 2_pInt**3_pInt, &
   CPFEM_RESTOREJACOBIAN = 2_pInt**4_pInt


 public :: &
   CPFEM_general, &
   CPFEM_initAll

contains


!--------------------------------------------------------------------------------------------------
!> @brief call (thread safe) all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll(el,ip)
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
#ifdef FEM
 use FEZoo, only: &
   FEZoo_init
#endif

 implicit none
 integer(pInt), intent(in) ::                        el, &                                          !< FE el number
                                                     ip                                             !< FE integration point number

 !$OMP CRITICAL (init)
   if (.not. CPFEM_init_done) then
#if defined(Spectral) || defined(FEM)
     call DAMASK_interface_init                                                                    ! Spectral and FEM interface to commandline
#endif
     call prec_init
     call IO_init
#ifdef FEM
     call FEZoo_init
#endif
     call numerics_init
     call debug_init
     call math_init
     call FE_init
     call mesh_init(ip, el)                                                                        ! pass on coordinates to alter calcMode of first ip
     call lattice_init
     call material_init
     call constitutive_init
     call crystallite_init
     call homogenization_init
     call CPFEM_init
#if defined(Marc4DAMASK) || defined(Abaqus)
     call DAMASK_interface_init                                                                    ! Spectral solver and FEM init is already done
#endif
     CPFEM_init_done = .true.
   endif
 !$OMP END CRITICAL (init)

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init
 use, intrinsic :: iso_fortran_env                                                                 ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use prec, only: &
   pInt
 use IO, only: &
   IO_read_realFile,&
   IO_read_intFile, &
   IO_timeStamp, &
   IO_error
 use numerics, only: &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive
 use FEsolving, only: &
#if defined(Marc4DAMASK) || defined(Abaqus)
   symmetricSolver, &
#endif
   restartRead, &
   modelName
#if defined(Marc4DAMASK) || defined(Abaqus)
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
#endif
 use material, only: &
   material_phase, &
   homogState, &
   phase_plasticity, &
   plasticState, &
   material_Nhomogenization
 use crystallite, only: &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Lp0, &
   crystallite_Fi0, &
   crystallite_Li0, &
   crystallite_dPdF0, &
   crystallite_Tstar0_v

 implicit none
 integer(pInt) :: k,l,m,ph,homog
 character(len=1024) :: rankStr

 mainProcess: if (worldrank == 0) then 
   write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
   write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 endif mainProcess

#if defined(Marc4DAMASK) || defined(Abaqus)
 ! initialize stress and jacobian to zero
 allocate(CPFEM_cs(6,mesh_maxNips,mesh_NcpElems)) ;                CPFEM_cs              = 0.0_pReal
 allocate(CPFEM_dcsdE(6,6,mesh_maxNips,mesh_NcpElems)) ;           CPFEM_dcsdE           = 0.0_pReal
 allocate(CPFEM_dcsdE_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dcsdE_knownGood = 0.0_pReal
#endif

 ! *** restore the last converged values of each essential variable from the binary file
 if (restartRead) then
   if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a)') '<< CPFEM >> restored state variables of last converged step from binary files'
     flush(6)
   endif
   
   write(rankStr,'(a1,i0)')'_',worldrank

   call IO_read_intFile(777,'recordedPhase'//trim(rankStr),modelName,size(material_phase))
   read (777,rec=1) material_phase
   close (777)

   call IO_read_realFile(777,'convergedF'//trim(rankStr),modelName,size(crystallite_F0))
   read (777,rec=1) crystallite_F0
   close (777)

   call IO_read_realFile(777,'convergedFp'//trim(rankStr),modelName,size(crystallite_Fp0))
   read (777,rec=1) crystallite_Fp0
   close (777)

   call IO_read_realFile(777,'convergedFi'//trim(rankStr),modelName,size(crystallite_Fi0))
   read (777,rec=1) crystallite_Fi0
   close (777)

   call IO_read_realFile(777,'convergedLp'//trim(rankStr),modelName,size(crystallite_Lp0))
   read (777,rec=1) crystallite_Lp0
   close (777)

   call IO_read_realFile(777,'convergedLi'//trim(rankStr),modelName,size(crystallite_Li0))
   read (777,rec=1) crystallite_Li0
   close (777)

   call IO_read_realFile(777,'convergeddPdF'//trim(rankStr),modelName,size(crystallite_dPdF0))
   read (777,rec=1) crystallite_dPdF0
   close (777)

   call IO_read_realFile(777,'convergedTstar'//trim(rankStr),modelName,size(crystallite_Tstar0_v))
   read (777,rec=1) crystallite_Tstar0_v
   close (777)

   call IO_read_realFile(777,'convergedStateConst'//trim(rankStr),modelName)
   m = 0_pInt
   readPlasticityInstances: do ph = 1_pInt, size(phase_plasticity)
     do k = 1_pInt, plasticState(ph)%sizeState
       do l = 1, size(plasticState(ph)%state0(1,:))
         m = m+1_pInt
         read(777,rec=m) plasticState(ph)%state0(k,l)
     enddo; enddo
   enddo readPlasticityInstances
   close (777)

   call IO_read_realFile(777,'convergedStateHomog'//trim(rankStr),modelName)
   m = 0_pInt
   readHomogInstances: do homog = 1_pInt, material_Nhomogenization
     do k = 1_pInt, homogState(homog)%sizeState
       do l = 1, size(homogState(homog)%state0(1,:))
         m = m+1_pInt
         read(777,rec=m) homogState(homog)%state0(k,l)
     enddo; enddo
   enddo readHomogInstances
   close (777)

#if defined(Marc4DAMASK) || defined(Abaqus)
   call IO_read_realFile(777,'convergeddcsdE',modelName,size(CPFEM_dcsdE))
   read (777,rec=1) CPFEM_dcsdE
   close (777)
#endif
   restartRead = .false.
 endif
#if defined(Marc4DAMASK) || defined(Abaqus)
 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0) then
   write(6,'(a32,1x,6(i8,1x))')   'CPFEM_cs:              ', shape(CPFEM_cs)
   write(6,'(a32,1x,6(i8,1x))')   'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
   write(6,'(a32,1x,6(i8,1x),/)') 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
   write(6,'(a32,l1)')            'symmetricSolver:       ', symmetricSolver
 endif
#endif
 flush(6)

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief perform initialization at first call, update variables and call the actual material model
!--------------------------------------------------------------------------------------------------
#if defined(Marc4DAMASK) || defined(Abaqus)
subroutine CPFEM_general(mode, parallelExecution, ffn, ffn1, temperature_inp, dt, elFE, ip, cauchyStress, jacobian)
#else
subroutine CPFEM_general(mode, ffn, ffn1, dt, elFE, ip)
#endif
 use numerics, only: &
   defgradTolerance, &
   iJacoStiffness, &
   worldrank
 use debug, only: &
   debug_level, &
   debug_CPFEM, &
   debug_levelBasic, &
   debug_levelExtensive, &
   debug_levelSelective, &
#if defined(Marc4DAMASK) || defined(Abaqus)
   debug_stressMaxLocation, &
   debug_stressMinLocation, &
   debug_jacobianMaxLocation, &
   debug_jacobianMinLocation, &
   debug_stressMax, &
   debug_stressMin, &
   debug_jacobianMax, &
   debug_jacobianMin, &
#endif
   debug_e, &
   debug_i
 use FEsolving, only: &
   terminallyIll, &
   FEsolving_execElem, &
   FEsolving_execIP, &
   restartWrite
 use math, only: &
   math_identity2nd, &
   math_mul33x33, &
   math_det33, &
   math_transpose33, &
   math_I3, &
   math_Mandel3333to66, &
   math_Mandel66to3333, &
   math_Mandel33to6, &
   math_Mandel6to33
 use mesh, only: &
   mesh_FEasCP, &
   mesh_NcpElems, &
   mesh_maxNips, &
   mesh_element
 use material, only: &
   microstructure_elemhomo, &
   plasticState, &
   sourceState, &
   homogState, &
   thermalState, &
   damageState, &
   vacancyfluxState, &
   hydrogenfluxState, &
   phaseAt, phasememberAt, &
   material_phase, &
   phase_plasticity, &
   temperature, &
   thermalMapping, &
   phase_Nsources, &
   material_homog, &
   material_Nhomogenization
 use crystallite, only: &
   crystallite_partionedF,&
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Fp, &
   crystallite_Fi0, &
   crystallite_Fi, &
   crystallite_Lp0, &
   crystallite_Lp, &
   crystallite_Li0, &
   crystallite_Li, &
   crystallite_dPdF0, &
   crystallite_dPdF, &
   crystallite_Tstar0_v, &
   crystallite_Tstar_v
 use homogenization, only: &  
   materialpoint_F, &
   materialpoint_F0, &
#if defined(Marc4DAMASK) || defined(Abaqus)
   materialpoint_P, &
   materialpoint_dPdF, &
   materialpoint_results, &
   materialpoint_sizeResults, &
#endif
   materialpoint_stressAndItsTangent, &
   materialpoint_postResults
 use IO, only: &
   IO_write_jobRealFile, &
   IO_warning
 use DAMASK_interface

 implicit none
 integer(pInt), intent(in) ::                        elFE, &                                        !< FE element number
                                                     ip                                             !< integration point number
 real(pReal), intent(in) ::                          dt                                             !< time increment
 real(pReal), dimension (3,3), intent(in) ::         ffn, &                                         !< deformation gradient for t=t0
                                                     ffn1                                           !< deformation gradient for t=t1
 integer(pInt), intent(in) ::                        mode                                           !< computation mode  1: regular computation plus aging of results
#if defined(Marc4DAMASK) || defined(Abaqus)
 real(pReal), intent(in) ::                          temperature_inp                                !< temperature
 logical, intent(in) ::                              parallelExecution                              !< flag indicating parallel computation of requested IPs
 real(pReal), dimension(6), intent(out) ::           cauchyStress                                   !< stress vector in Mandel notation
 real(pReal), dimension(6,6), intent(out) ::         jacobian                                       !< jacobian in Mandel notation (Consistent tangent dcs/dE)

 real(pReal)                                         J_inverse, &                                   ! inverse of Jacobian
                                                     rnd
 real(pReal), dimension (3,3) ::                     Kirchhoff, &                                   ! Piola-Kirchhoff stress in Matrix notation
                                                     cauchyStress33                                 ! stress vector in Matrix notation
 real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                     H, &
                                                     jacobian3333                                   ! jacobian in Matrix notation
#else
 logical, parameter ::                               parallelExecution = .true.
#endif

 integer(pInt)                                       elCP, &                                        ! crystal plasticity element number
                                                     i, j, k, l, m, n, ph, homog, mySource
 logical                                             updateJaco                                     ! flag indicating if JAcobian has to be updated
 character(len=1024) :: rankStr

#if defined(Marc4DAMASK) || defined(Abaqus)
 elCP = mesh_FEasCP('elem',elFE)
#else
 elCP = elFE
#endif

 if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt &
     .and. elCP == debug_e .and. ip == debug_i) then
   write(6,'(/,a)') '#############################################'
   write(6,'(a1,a22,1x,i8,a13)')   '#','element',        elCP,         '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','ip',             ip,           '#'
   write(6,'(a1,a22,1x,f15.7,a6)') '#','theTime',        theTime,      '#'
   write(6,'(a1,a22,1x,f15.7,a6)') '#','theDelta',       theDelta,     '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','theInc',         theInc,       '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','cycleCounter',   cycleCounter, '#'
   write(6,'(a1,a22,1x,i8,a13)')   '#','computationMode',mode,         '#'
   if (terminallyIll) &
   write(6,'(a,/)') '#           --- terminallyIll ---           #'
   write(6,'(a,/)') '#############################################'; flush (6)
 endif


#if defined(Marc4DAMASK) || defined(Abaqus)
 if (iand(mode, CPFEM_BACKUPJACOBIAN) /= 0_pInt) &
   CPFEM_dcsde_knownGood = CPFEM_dcsde
 if (iand(mode, CPFEM_RESTOREJACOBIAN) /= 0_pInt) &
   CPFEM_dcsde = CPFEM_dcsde_knownGood
#endif

 !*** age results and write restart data if requested
 if (iand(mode, CPFEM_AGERESULTS) /= 0_pInt) then
   crystallite_F0  = crystallite_partionedF                                                    ! crystallite deformation (_subF is perturbed...)
   crystallite_Fp0 = crystallite_Fp                                                            ! crystallite plastic deformation
   crystallite_Lp0 = crystallite_Lp                                                            ! crystallite plastic velocity
   crystallite_Fi0 = crystallite_Fi                                                            ! crystallite intermediate deformation
   crystallite_Li0 = crystallite_Li                                                            ! crystallite intermediate velocity
   crystallite_dPdF0 = crystallite_dPdF                                                        ! crystallite stiffness
   crystallite_Tstar0_v = crystallite_Tstar_v                                                  ! crystallite 2nd Piola Kirchhoff stress

   forall ( i = 1:size(plasticState    )) plasticState(i)%state0     = plasticState(i)%state   ! copy state in this lenghty way because: A component cannot be an array if the encompassing structure is an array
   do i = 1, size(sourceState)
     do mySource = 1,phase_Nsources(i)
       sourceState(i)%p(mySource)%state0 = sourceState(i)%p(mySource)%state                    ! copy state in this lenghty way because: A component cannot be an array if the encompassing structure is an array
   enddo; enddo
   if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) then
     write(6,'(a)') '<< CPFEM >> aging states'
     if (debug_e <= mesh_NcpElems .and. debug_i <= mesh_maxNips) then
       write(6,'(a,1x,i8,1x,i2,1x,i4,/,(12x,6(e20.8,1x)),/)') &
             '<< CPFEM >> aged state of elFE ip grain',debug_e, debug_i, 1, &
              plasticState(phaseAt(1,debug_i,debug_e))%state(:,phasememberAt(1,debug_i,debug_e))
       endif
   endif

   do homog = 1_pInt, material_Nhomogenization
     homogState       (homog)%state0 =  homogState       (homog)%state
     thermalState     (homog)%state0 =  thermalState     (homog)%state
     damageState      (homog)%state0 =  damageState      (homog)%state
     vacancyfluxState (homog)%state0 =  vacancyfluxState (homog)%state
     hydrogenfluxState(homog)%state0 =  hydrogenfluxState(homog)%state
   enddo
  

   ! * dump the last converged values of each essential variable to a binary file

   if (restartWrite) then
     if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt) &
        write(6,'(a)') '<< CPFEM >> writing state variables of last converged step to binary files'
   
     write(rankStr,'(a1,i0)')'_',worldrank

     call IO_write_jobRealFile(777,'recordedPhase'//trim(rankStr),size(material_phase))
     write (777,rec=1) material_phase
     close (777)

     call IO_write_jobRealFile(777,'convergedF'//trim(rankStr),size(crystallite_F0))
     write (777,rec=1) crystallite_F0
     close (777)

     call IO_write_jobRealFile(777,'convergedFp'//trim(rankStr),size(crystallite_Fp0))
     write (777,rec=1) crystallite_Fp0
     close (777)

     call IO_write_jobRealFile(777,'convergedFi'//trim(rankStr),size(crystallite_Fi0))
     write (777,rec=1) crystallite_Fi0
     close (777)

     call IO_write_jobRealFile(777,'convergedLp'//trim(rankStr),size(crystallite_Lp0))
     write (777,rec=1) crystallite_Lp0
     close (777)

     call IO_write_jobRealFile(777,'convergedLi'//trim(rankStr),size(crystallite_Li0))
     write (777,rec=1) crystallite_Li0
     close (777)

     call IO_write_jobRealFile(777,'convergeddPdF'//trim(rankStr),size(crystallite_dPdF0))
     write (777,rec=1) crystallite_dPdF0
     close (777)

     call IO_write_jobRealFile(777,'convergedTstar'//trim(rankStr),size(crystallite_Tstar0_v))
     write (777,rec=1) crystallite_Tstar0_v
     close (777)

     call IO_write_jobRealFile(777,'convergedStateConst'//trim(rankStr))
     m = 0_pInt
     writePlasticityInstances: do ph = 1_pInt, size(phase_plasticity)
       do k = 1_pInt, plasticState(ph)%sizeState
         do l = 1, size(plasticState(ph)%state0(1,:))
           m = m+1_pInt
           write(777,rec=m) plasticState(ph)%state0(k,l)
       enddo; enddo
     enddo writePlasticityInstances
     close (777)

     call IO_write_jobRealFile(777,'convergedStateHomog'//trim(rankStr))
     m = 0_pInt
     writeHomogInstances: do homog = 1_pInt, material_Nhomogenization
       do k = 1_pInt, homogState(homog)%sizeState
         do l = 1, size(homogState(homog)%state0(1,:))
           m = m+1_pInt
           write(777,rec=m) homogState(homog)%state0(k,l)
       enddo; enddo
     enddo writeHomogInstances
     close (777)

#if defined(Marc4DAMASK) || defined(Abaqus)
     call IO_write_jobRealFile(777,'convergeddcsdE',size(CPFEM_dcsdE))
     write (777,rec=1) CPFEM_dcsdE
     close (777)
#endif

   endif
 endif                                                              ! results aging



 !*** collection of FEM input with returning of randomize odd stress and jacobian
 !*   If no parallel execution is required, there is no need to collect FEM input

 if (.not. parallelExecution) then
#if defined(Marc4DAMASK) || defined(Abaqus)
   temperature(material_homog(ip,elCP))%p(thermalMapping(material_homog(ip,elCP))%p(ip,elCP)) = &
     temperature_inp
#endif
   materialpoint_F0(1:3,1:3,ip,elCP) = ffn
   materialpoint_F(1:3,1:3,ip,elCP) = ffn1

 elseif (iand(mode, CPFEM_COLLECT) /= 0_pInt) then
#if defined(Marc4DAMASK) || defined(Abaqus)
   call random_number(rnd)
   if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
   CPFEM_cs(1:6,ip,elCP) = rnd * CPFEM_odd_stress
   CPFEM_dcsde(1:6,1:6,ip,elCP) = CPFEM_odd_jacobian * math_identity2nd(6)
   temperature(material_homog(ip,elCP))%p(thermalMapping(material_homog(ip,elCP))%p(ip,elCP)) = &
     temperature_inp
#endif
   materialpoint_F0(1:3,1:3,ip,elCP) = ffn
   materialpoint_F(1:3,1:3,ip,elCP) = ffn1
   CPFEM_calc_done = .false.
 endif                                                              ! collection



 !*** calculation of stress and jacobian

 if (iand(mode, CPFEM_CALCRESULTS) /= 0_pInt) then

   !*** deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
   validCalculation: if (terminallyIll &
                    .or. outdatedFFN1 &
                    .or. any(abs(ffn1 - materialpoint_F(1:3,1:3,ip,elCP)) > defgradTolerance)) then
     if (any(abs(ffn1 - materialpoint_F(1:3,1:3,ip,elCP)) > defgradTolerance)) then
       if (iand(debug_level(debug_CPFEM), debug_levelBasic) /=  0_pInt) then
           write(6,'(a,1x,i8,1x,i2)') '<< CPFEM >> OUTDATED at elFE ip',elFE,ip
           write(6,'(a,/,3(12x,3(f10.6,1x),/))') '<< CPFEM >> FFN1 old:',&
                                             math_transpose33(materialpoint_F(1:3,1:3,ip,elCP))
           write(6,'(a,/,3(12x,3(f10.6,1x),/))') '<< CPFEM >> FFN1 now:',math_transpose33(ffn1)
       endif
       outdatedFFN1 = .true.
     endif
#if defined(Marc4DAMASK) || defined(Abaqus)
     call random_number(rnd)
     if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
     CPFEM_cs(1:6,ip,elCP) = rnd*CPFEM_odd_stress
     CPFEM_dcsde(1:6,1:6,ip,elCP) = CPFEM_odd_jacobian*math_identity2nd(6)
#endif

   !*** deformation gradient is not outdated

   else validCalculation
     updateJaco = mod(cycleCounter,iJacoStiffness) == 0
     !* no parallel computation, so we use just one single elFE and ip for computation

     if (.not. parallelExecution) then
       FEsolving_execElem(1)     = elCP
       FEsolving_execElem(2)     = elCP
       if (.not. microstructure_elemhomo(mesh_element(4,elCP)) .or. &                               ! calculate unless homogeneous
                (microstructure_elemhomo(mesh_element(4,elCP)) .and. ip == 1_pInt)) then            ! and then only first ip
         FEsolving_execIP(1,elCP) = ip
         FEsolving_execIP(2,elCP) = ip
         if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /=  0_pInt) &
           write(6,'(a,i8,1x,i2)') '<< CPFEM >> calculation for elFE ip ',elFE,ip
         call materialpoint_stressAndItsTangent(updateJaco, dt)                                     ! calculate stress and its tangent
         call materialpoint_postResults()
       endif

     !* parallel computation and calulation not yet done

     elseif (.not. CPFEM_calc_done) then
       if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) &
         write(6,'(a,i8,a,i8)') '<< CPFEM >> calculation for elements ',FEsolving_execElem(1),&
                                                                 ' to ',FEsolving_execElem(2)
       call materialpoint_stressAndItsTangent(updateJaco, dt)                                       ! calculate stress and its tangent (parallel execution inside)
       call materialpoint_postResults()
       CPFEM_calc_done = .true.
     endif

     !* map stress and stiffness (or return odd values if terminally ill)
#if defined(Marc4DAMASK) || defined(Abaqus)
     terminalIllness: if ( terminallyIll ) then

       call random_number(rnd)
       if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
       CPFEM_cs(1:6,ip,elCP) = rnd * CPFEM_odd_stress
       CPFEM_dcsde(1:6,1:6,ip,elCP) = CPFEM_odd_jacobian * math_identity2nd(6)

     else terminalIllness

       if (microstructure_elemhomo(mesh_element(4,elCP)) .and. ip > 1_pInt) then                    ! me homogenous? --> copy from first ip
         materialpoint_P(1:3,1:3,ip,elCP) = materialpoint_P(1:3,1:3,1,elCP)
         materialpoint_F(1:3,1:3,ip,elCP) = materialpoint_F(1:3,1:3,1,elCP)
         materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,elCP) = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,elCP)
         materialpoint_results(1:materialpoint_sizeResults,ip,elCP) = &
         materialpoint_results(1:materialpoint_sizeResults,1,elCP)
       endif

       ! translate from P to CS
       Kirchhoff = math_mul33x33(materialpoint_P(1:3,1:3,ip,elCP), math_transpose33(materialpoint_F(1:3,1:3,ip,elCP)))
       J_inverse  = 1.0_pReal / math_det33(materialpoint_F(1:3,1:3,ip,elCP))
       CPFEM_cs(1:6,ip,elCP) = math_Mandel33to6(J_inverse * Kirchhoff)

       !  translate from dP/dF to dCS/dE
       H = 0.0_pReal
       do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
         H(i,j,k,l) = H(i,j,k,l) + &
                       materialpoint_F(j,m,ip,elCP) * &
                       materialpoint_F(l,n,ip,elCP) * &
                       materialpoint_dPdF(i,m,k,n,ip,elCP) - &
                       math_I3(j,l) * materialpoint_F(i,m,ip,elCP) * materialpoint_P(k,m,ip,elCP) + &
                       0.5_pReal * (math_I3(i,k) * Kirchhoff(j,l) + math_I3(j,l) * Kirchhoff(i,k) + &
                                  math_I3(i,l) * Kirchhoff(j,k) + math_I3(j,k) * Kirchhoff(i,l))
       enddo; enddo; enddo; enddo; enddo; enddo

       forall(i=1:3, j=1:3,k=1:3,l=1:3) &
         H_sym(i,j,k,l) = 0.25_pReal * (H(i,j,k,l) + H(j,i,k,l) + H(i,j,l,k) + H(j,i,l,k))

       CPFEM_dcsde(1:6,1:6,ip,elCP) = math_Mandel3333to66(J_inverse * H_sym)

     endif terminalIllness
#endif

   endif validCalculation

#if defined(Marc4DAMASK) || defined(Abaqus)
   !* report stress and stiffness
   if ((iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) &
        .and. ((debug_e == elCP .and. debug_i == ip) &
               .or. .not. iand(debug_level(debug_CPFEM), debug_levelSelective) /= 0_pInt)) then
       write(6,'(a,i8,1x,i2,/,12x,6(f10.3,1x)/)') &
         '<< CPFEM >> stress/MPa at elFE ip ',   elFE, ip, CPFEM_cs(1:6,ip,elCP)*1.0e-6_pReal
       write(6,'(a,i8,1x,i2,/,6(12x,6(f10.3,1x)/))') &
         '<< CPFEM >> Jacobian/GPa at elFE ip ', elFE, ip, transpose(CPFEM_dcsdE(1:6,1:6,ip,elCP))*1.0e-9_pReal
       flush(6)
   endif
#endif

 endif

#if defined(Marc4DAMASK) || defined(Abaqus)
 !*** warn if stiffness close to zero
 if (all(abs(CPFEM_dcsdE(1:6,1:6,ip,elCP)) < 1e-10_pReal)) call IO_warning(601,elCP,ip)

 !*** copy to output if using commercial FEM solver
 cauchyStress = CPFEM_cs   (1:6,    ip,elCP)
 jacobian     = CPFEM_dcsdE(1:6,1:6,ip,elCP)


 !*** remember extreme values of stress ...
 cauchyStress33 = math_Mandel6to33(CPFEM_cs(1:6,ip,elCP))
 if (maxval(cauchyStress33) > debug_stressMax) then
   debug_stressMaxLocation = [elCP, ip]
   debug_stressMax = maxval(cauchyStress33)
 endif
 if (minval(cauchyStress33) < debug_stressMin) then
   debug_stressMinLocation = [elCP, ip]
   debug_stressMin = minval(cauchyStress33)
 endif
 !*** ... and Jacobian
 jacobian3333 = math_Mandel66to3333(CPFEM_dcsdE(1:6,1:6,ip,elCP))
 if (maxval(jacobian3333) > debug_jacobianMax) then
   debug_jacobianMaxLocation = [elCP, ip]
   debug_jacobianMax = maxval(jacobian3333)
 endif
 if (minval(jacobian3333) < debug_jacobianMin) then
   debug_jacobianMinLocation = [elCP, ip]
   debug_jacobianMin = minval(jacobian3333)
 endif
#endif

end subroutine CPFEM_general

end module CPFEM
