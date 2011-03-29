!* $Id$
!##############################################################
MODULE CPFEM
!##############################################################
!    *** CPFEM engine ***
!
use prec, only:                                   pReal, &
                                                  pInt
implicit none
 
real(pReal), parameter ::                         CPFEM_odd_stress    = 1e15_pReal, &
                                                  CPFEM_odd_jacobian  = 1e50_pReal

real(pReal), dimension (:,:,:),   allocatable ::  CPFEM_cs                            ! Cauchy stress
real(pReal), dimension (:,:,:,:), allocatable ::  CPFEM_dcsdE                         ! Cauchy stress tangent
real(pReal), dimension (:,:,:,:), allocatable ::  CPFEM_dcsdE_knownGood               ! known good tangent

logical ::                                        CPFEM_init_done       = .false., &  ! remember whether init has been done already
                                                  CPFEM_init_inProgress = .false., &  ! remember whether first IP is currently performing init
                                                  CPFEM_calc_done       = .false.     ! remember whether first IP has already calced the results
 

CONTAINS

!*********************************************************
!***    call (thread safe) all module initializations  ***
!*********************************************************

subroutine CPFEM_initAll(Temperature,element,IP)

  use prec, only:                                     pReal, &
                                                      prec_init
  use numerics, only:                                 numerics_init
  use debug, only:                                    debug_init
  use FEsolving, only:                                FE_init
  use math, only:                                     math_init
  use mesh, only:                                     mesh_init
  use lattice, only:                                  lattice_init
  use material, only:                                 material_init
  use constitutive, only:                             constitutive_init
  use crystallite, only:                              crystallite_init
  use homogenization, only:                           homogenization_init
  use IO, only:                                       IO_init
  use mpie_interface
  implicit none

  integer(pInt), intent(in) ::                        element, &          ! FE element number
                                                      IP                  ! FE integration point number
  real(pReal), intent(in) ::                          Temperature         ! temperature
  real(pReal)   rnd
  integer(pInt) i,n

  ! initialization step (three dimensional stress state check missing?)

  if (.not. CPFEM_init_done) then
    call random_number(rnd)
    do i=1,int(256.0*rnd)
      n = n+1_pInt                                                      ! wasting random amount of time...
    enddo                                                               ! ...to break potential race in multithreading
    n = n+1_pInt
    if (.not. CPFEM_init_inProgress) then                               ! yes my thread won!
      CPFEM_init_inProgress = .true.
      call prec_init()
      call IO_init()
      call numerics_init()
      call debug_init()
      call math_init()
      call FE_init()
      call mesh_init(IP, element)                ! pass on coordinates to alter calcMode of first ip
      call lattice_init()
      call material_init()
      call constitutive_init()
      call crystallite_init(Temperature)         ! (have to) use temperature of first IP for whole model
      call homogenization_init(Temperature)
      call CPFEM_init()
      call mpie_interface_init()
      CPFEM_init_done = .true.
      CPFEM_init_inProgress = .false.
    else                                                                ! loser, loser...
      do while (CPFEM_init_inProgress)
      enddo
    endif
  endif

end subroutine

!*********************************************************
!***    allocate the arrays defined in module CPFEM    ***
!***    and initialize them                            ***
!*********************************************************

subroutine CPFEM_init()

  use prec, only:                                 pInt
  use debug, only:                                debug_verbosity
  use IO, only:                                   IO_read_jobBinaryFile
  use FEsolving, only:                            parallelExecution, &
                                                  symmetricSolver, &
                                                  restartRead, &
                                                  restartJob
  use mesh, only:                                 mesh_NcpElems, &
                                                  mesh_maxNips
  use material, only:                             homogenization_maxNgrains, &
                                                  material_phase
  use constitutive, only:                         constitutive_state0
  use crystallite, only:                          crystallite_F0, &
                                                  crystallite_Fp0, &
                                                  crystallite_Lp0, &
                                                  crystallite_dPdF0, &
                                                  crystallite_Tstar0_v
  use homogenization, only:                       homogenization_sizeState, &
                                                  homogenization_state0, &
                                                  materialpoint_F, &
                                                  materialpoint_F0


  implicit none

  integer(pInt) i,j,k,l,m

  ! initialize stress and jacobian to zero 
  allocate(CPFEM_cs(6,mesh_maxNips,mesh_NcpElems)) ;                CPFEM_cs              = 0.0_pReal
  allocate(CPFEM_dcsdE(6,6,mesh_maxNips,mesh_NcpElems)) ;           CPFEM_dcsdE           = 0.0_pReal
  allocate(CPFEM_dcsdE_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dcsdE_knownGood = 0.0_pReal

  ! *** restore the last converged values of each essential variable from the binary file
  if (restartRead) then
    if (debug_verbosity > 0) then
      !$OMP CRITICAL (write2out)
       write(6,'(a)') '<< CPFEM >> Restored state variables of last converged step from binary files'
      !$OMP END CRITICAL (write2out)
    endif
    if (IO_read_jobBinaryFile(777,'recordedPhase',restartJob,size(material_phase))) then
      read (777,rec=1) material_phase
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergedF',restartJob,size(crystallite_F0))) then
      read (777,rec=1) crystallite_F0
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergedFp',restartJob,size(crystallite_Fp0))) then
      read (777,rec=1) crystallite_Fp0
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergedLp',restartJob,size(crystallite_Lp0))) then
      read (777,rec=1) crystallite_Lp0
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergeddPdF',restartJob,size(crystallite_dPdF0))) then
      read (777,rec=1) crystallite_dPdF0
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergedTstar',restartJob,size(crystallite_Tstar0_v))) then
      read (777,rec=1) crystallite_Tstar0_v
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergedStateConst',restartJob)) then
      m = 0_pInt
      do i = 1,homogenization_maxNgrains; do j = 1,mesh_maxNips; do k = 1,mesh_NcpElems
        do l = 1,size(constitutive_state0(i,j,k)%p)
          m = m+1_pInt
          read(777,rec=m) constitutive_state0(i,j,k)%p(l)
        enddo
      enddo; enddo; enddo
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergedStateHomog',restartJob)) then
      m = 0_pInt
      do k = 1,mesh_NcpElems; do j = 1,mesh_maxNips
        do l = 1,homogenization_sizeState(j,k)
          m = m+1_pInt
          read(777,rec=m) homogenization_state0(j,k)%p(l)
        enddo
      enddo; enddo
      close (777)
    endif
    if (IO_read_jobBinaryFile(777,'convergeddcsdE',restartJob,size(CPFEM_dcsdE))) then
      read (777,rec=1) CPFEM_dcsdE
      close (777)
    endif
    restartRead = .false.
  endif
  ! *** end of restoring

  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,*) '<<<+-  cpfem init  -+>>>'
    write(6,*) '$Id$'
    write(6,*)
    if (debug_verbosity > 0) then
      write(6,'(a32,x,6(i5,x))') 'CPFEM_cs:              ', shape(CPFEM_cs)
      write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
      write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
      write(6,*)
      write(6,*) 'parallelExecution:    ', parallelExecution
      write(6,*) 'symmetricSolver:      ', symmetricSolver
    endif
    call flush(6)
  !$OMP END CRITICAL (write2out)

endsubroutine


!***********************************************************************
!***    perform initialization at first call, update variables and   ***
!***    call the actual material model                               ***
!***********************************************************************
subroutine CPFEM_general(mode, ffn, ffn1, Temperature, dt, element, IP, cauchyStress,&
      & jacobian, pstress, dPdF)
  ! note: cauchyStress = Cauchy stress cs(6) and jacobian = Consistent tangent dcs/dE

  !*** variables and functions from other modules ***!
  use prec, only:                                     pReal, &
                                                      pInt
  use numerics, only:                                 relevantStrain, &
                                                      defgradTolerance, &
                                                      iJacoStiffness
  use debug, only:                                    debug_e, &
                                                      debug_i, &
                                                      debug_g, &
                                                      debug_selectiveDebugger, &
                                                      debug_verbosity, &
                                                      debug_stressMaxLocation, &
                                                      debug_stressMinLocation, &
                                                      debug_jacobianMaxLocation, &
                                                      debug_jacobianMinLocation, &
                                                      debug_stressMax, &
                                                      debug_stressMin, &
                                                      debug_jacobianMax, &
                                                      debug_jacobianMin
  use FEsolving, only:                                parallelExecution, &
                                                      outdatedFFN1, &
                                                      terminallyIll, &
                                                      cycleCounter, &
                                                      theInc, &
                                                      theTime, &
                                                      theDelta, &
                                                      FEsolving_execElem, &
                                                      FEsolving_execIP, &
                                                      restartWrite
  use math, only:                                     math_identity2nd, &
                                                      math_mul33x33, &
                                                      math_det3x3, &
                                                      math_transpose3x3, &
                                                      math_I3, &
                                                      math_Mandel3333to66, &
                                                      math_Mandel66to3333, &
                                                      math_Mandel33to6, &
                                                      math_Mandel6to33
  use mesh, only:                                     mesh_FEasCP, &
                                                      mesh_NcpElems, &
                                                      mesh_maxNips, &
                                                      mesh_element, &
                                                      FE_Nips
  use material, only:                                 homogenization_maxNgrains, &
                                                      microstructure_elemhomo, &
                                                      material_phase
  use constitutive, only:                             constitutive_state0,constitutive_state
  use crystallite, only:                              crystallite_F0, &
                                                      crystallite_partionedF, &
                                                      crystallite_Fp0, &
                                                      crystallite_Fp, &
                                                      crystallite_Lp0, &
                                                      crystallite_Lp, &
                                                      crystallite_dPdF0, &
                                                      crystallite_dPdF, &
                                                      crystallite_Tstar0_v, &
                                                      crystallite_Tstar_v
  use homogenization, only:                           homogenization_sizeState, &
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
  use IO, only:                                       IO_write_jobBinaryFile, &
                                                      IO_warning
  use mpie_interface
  
  implicit none
  
  !*** input variables ***!
  integer(pInt), intent(in) ::                        element, &          ! FE element number
                                                      IP                  ! FE integration point number
  real(pReal), intent(inout) ::                       Temperature         ! temperature
  real(pReal), intent(in) ::                          dt                  ! time increment
  real(pReal), dimension (3,3), intent(in) ::         ffn, &              ! deformation gradient for t=t0
                                                      ffn1                ! deformation gradient for t=t1
  integer(pInt), intent(in) ::                        mode                ! computation mode  1: regular computation plus aging of results
                                                                          !                   2: regular computation
                                                                          !                   3: collection of FEM data
                                                                          !                   4: backup tangent from former converged inc
                                                                          !                   5: restore tangent from former converged inc
                                                                          !                   6: recycling of former results (MARC speciality)
  
  !*** output variables ***!
  real(pReal), dimension(6), intent(out) ::           cauchyStress        ! stress vector in Mandel notation
  real(pReal), dimension(6,6), intent(out) ::         jacobian            ! jacobian in Mandel notation
  real(pReal), dimension (3,3), intent(out) ::        pstress             ! Piola-Kirchhoff stress in Matrix notation
  real(pReal), dimension (3,3,3,3), intent(out) ::    dPdF                ! 
                 
  !*** local variables ***!
  real(pReal)                                         J_inverse, &        ! inverse of Jacobian
                                                      rnd
  real(pReal), dimension (3,3) ::                     Kirchhoff, &        ! Piola-Kirchhoff stress in Matrix notation
                                                      cauchyStress33      ! stress vector in Matrix notation
  real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                      H, &
                                                      jacobian3333        ! jacobian in Matrix notation
  integer(pInt)                                       cp_en, &            ! crystal plasticity element number
                                                      i, &
                                                      j, &
                                                      k, &
                                                      l, &
                                                      m, &
                                                      n, &
                                                      e
  logical                                             updateJaco          ! flag indicating if JAcobian has to be updated
  
  !*** global variables ***!
  ! CPFEM_cs, &
  ! CPFEM_dcsdE, &
  ! CPFEM_dcsdE_knownGood, &
  ! CPFEM_init_done, &
  ! CPFEM_calc_done, &
  ! CPFEM_odd_stress, &
  ! CPFEM_odd_jacobian
  
  cp_en = mesh_FEasCP('elem',element)
  
  if (debug_verbosity > 0 .and. cp_en == 1 .and. IP == 1) then
    !$OMP CRITICAL (write2out)
      write(6,*)
      write(6,'(a)') '#############################################'
      write(6,'(a1,a22,x,f15.7,a6)') '#','theTime',theTime,'#'
      write(6,'(a1,a22,x,f15.7,a6)') '#','theDelta',theDelta,'#'
      write(6,'(a1,a22,x,i8,a13)') '#','theInc',theInc,'#'
      write(6,'(a1,a22,x,i8,a13)') '#','cycleCounter',cycleCounter,'#'
      write(6,'(a1,a22,x,i8,a13)') '#','computationMode',mode,'#'
      write(6,'(a)') '#############################################'
      write(6,*)
      call flush (6)
    !$OMP END CRITICAL (write2out)
  endif

  ! according to our "mode" we decide what to do
  select case (mode)
    
    ! --+>> REGULAR COMPUTATION (WITH AGING OF RESULTS IF MODE == 1) <<+-- 
    case (1,2,8,9)
      ! age results if mode == 1
      if (mode == 1 .or. mode == 8) then
        crystallite_F0  = crystallite_partionedF                          ! crystallite deformation (_subF is perturbed...)
        crystallite_Fp0 = crystallite_Fp                                  ! crystallite plastic deformation
        crystallite_Lp0 = crystallite_Lp                                  ! crystallite plastic velocity
        crystallite_dPdF0 = crystallite_dPdF                              ! crystallite stiffness
        crystallite_Tstar0_v = crystallite_Tstar_v                        ! crystallite 2nd Piola Kirchhoff stress 
        forall ( i = 1:homogenization_maxNgrains, &
                 j = 1:mesh_maxNips, &
                 k = 1:mesh_NcpElems ) &
          constitutive_state0(i,j,k)%p = constitutive_state(i,j,k)%p      ! microstructure of crystallites
        if (debug_verbosity > 0) then
          !$OMP CRITICAL (write2out)
            write(6,'(a)') '<< CPFEM >> Aging states'
            if (debug_e == cp_en .and. debug_i == IP) then
              write(6,'(a,x,i5,x,i2,x,i3,/,(12(x),6(e20.8,x),/))') '<< CPFEM >> AGED state of element ip grain',&
                                                              cp_en, IP, 1, constitutive_state(1,IP,cp_en)%p
            endif
          !$OMP END CRITICAL (write2out)
        endif
        !$OMP PARALLEL DO
          do k = 1,mesh_NcpElems
            do j = 1,mesh_maxNips
              if (homogenization_sizeState(j,k) > 0_pInt) &
                homogenization_state0(j,k)%p = homogenization_state(j,k)%p  ! internal state of homogenization scheme
            enddo
          enddo
        !$OMP END PARALLEL DO


        ! *** dump the last converged values of each essential variable to a binary file
        if (restartWrite) then 
          if (debug_verbosity > 0) then
           !$OMP CRITICAL (write2out)
             write(6,'(a)') '<< CPFEM >> Writing state variables of last converged step to binary files'
           !$OMP END CRITICAL (write2out)
          endif
          if (IO_write_jobBinaryFile(777,'recordedPhase',size(material_phase))) then
            write (777,rec=1) material_phase
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergedF',size(crystallite_F0))) then
            write (777,rec=1) crystallite_F0
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergedFp',size(crystallite_Fp0))) then
            write (777,rec=1) crystallite_Fp0
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergedLp',size(crystallite_Lp0))) then
            write (777,rec=1) crystallite_Lp0
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergeddPdF',size(crystallite_dPdF0))) then
            write (777,rec=1) crystallite_dPdF0
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergedTstar',size(crystallite_Tstar0_v))) then
            write (777,rec=1) crystallite_Tstar0_v
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergedStateConst')) then
            m = 0_pInt
            do i = 1,homogenization_maxNgrains; do j = 1,mesh_maxNips; do k = 1,mesh_NcpElems
              do l = 1,size(constitutive_state0(i,j,k)%p)
                m = m+1_pInt
                write(777,rec=m) constitutive_state0(i,j,k)%p(l)
              enddo
            enddo; enddo; enddo
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergedStateHomog')) then
            m = 0_pInt
            do k = 1,mesh_NcpElems; do j = 1,mesh_maxNips
              do l = 1,homogenization_sizeState(j,k)
                m = m+1_pInt
                write(777,rec=m) homogenization_state0(j,k)%p(l)
              enddo
            enddo; enddo
            close (777)
          endif
          if (IO_write_jobBinaryFile(777,'convergeddcsdE',size(CPFEM_dcsdE))) then
            write (777,rec=1) CPFEM_dcsdE
            close (777)
          endif
        endif
        ! *** end of dumping
      endif

      if (mode == 8 .or. mode == 9) then                                  ! Abaqus explicit skips collect
        materialpoint_Temperature(IP,cp_en) = Temperature
        materialpoint_F0(1:3,1:3,IP,cp_en) = ffn
        materialpoint_F(1:3,1:3,IP,cp_en) = ffn1
      endif

      ! deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
      if (terminallyIll .or. outdatedFFN1 .or. any(abs(ffn1 - materialpoint_F(1:3,1:3,IP,cp_en)) > defgradTolerance)) then
        if (.not. terminallyIll .and. .not. outdatedFFN1) then 
          if (debug_verbosity > 0) then
            !$OMP CRITICAL (write2out)
              write(6,'(a,x,i5,x,i2)') '<< CPFEM >> OUTDATED at element ip',cp_en,IP
              write(6,'(a,/,3(12(x),3(f10.6,x),/))') '<< CPFEM >> FFN1 old:',math_transpose3x3(materialpoint_F(1:3,1:3,IP,cp_en))
              write(6,'(a,/,3(12(x),3(f10.6,x),/))') '<< CPFEM >> FFN1 now:',math_transpose3x3(ffn1)
            !$OMP END CRITICAL (write2out)
          endif
          outdatedFFN1 = .true.
        endif
        call random_number(rnd)
        rnd = 2.0_pReal * rnd - 1.0_pReal
        CPFEM_cs(1:6,IP,cp_en) = rnd*CPFEM_odd_stress
        CPFEM_dcsde(1:6,1:6,IP,cp_en) = CPFEM_odd_jacobian*math_identity2nd(6)
      
      ! deformation gradient is not outdated
      else
        ! set flag for Jacobian update
        updateJaco = mod(cycleCounter,iJacoStiffness) == 0
        
        ! no parallel computation
        if (.not. parallelExecution) then
          ! we just take one single element and IP
          FEsolving_execElem(1)     = cp_en
          FEsolving_execElem(2)     = cp_en
          FEsolving_execIP(1,cp_en) = IP
          FEsolving_execIP(2,cp_en) = IP
          if (debug_verbosity > 0) then
            !$OMP CRITICAL (write2out)
              write(6,'(a,i5,x,i2)') '<< CPFEM >> Calculation for element ip ',cp_en,IP
            !$OMP END CRITICAL (write2out)
          endif
          call materialpoint_stressAndItsTangent(updateJaco, dt)          ! calculate stress and its tangent
          call materialpoint_postResults(dt)                              ! post results
          
        ! parallel computation and calulation not yet done
        elseif (.not. CPFEM_calc_done) then
          if (debug_verbosity > 0) then
            !$OMP CRITICAL (write2out)
              write(6,'(a,i5,a,i5)') '<< CPFEM >> Calculation for elements ',FEsolving_execElem(1),' to ',FEsolving_execElem(2)
            !$OMP END CRITICAL (write2out)
          endif
          call materialpoint_stressAndItsTangent(updateJaco, dt)          ! calculate stress and its tangent (parallel execution inside)
          call materialpoint_postResults(dt)                              ! post results
          !$OMP PARALLEL DO
            do e = FEsolving_execElem(1),FEsolving_execElem(2)            ! loop over all parallely processed elements
              if (microstructure_elemhomo(mesh_element(4,e))) then        ! dealing with homogeneous element?
                forall (i = 2:FE_Nips(mesh_element(2,e)))                 ! copy results of first IP to all others
                  materialpoint_P(1:3,1:3,i,e) = materialpoint_P(1:3,1:3,1,e) 
                  materialpoint_F(1:3,1:3,i,e) = materialpoint_F(1:3,1:3,1,e) 
                  materialpoint_dPdF(1:3,1:3,1:3,1:3,i,e) = materialpoint_dPdF(1:3,1:3,1:3,1:3,1,e)
                  materialpoint_results(1:materialpoint_sizeResults,i,e) = materialpoint_results(1:materialpoint_sizeResults,1,e)
                end forall
              endif
            enddo
          !$OMP END PARALLEL DO
          CPFEM_calc_done = .true.
        endif
        
        if ( terminallyIll ) then
          call random_number(rnd)
          rnd = 2.0_pReal * rnd - 1.0_pReal
          CPFEM_cs(1:6,IP,cp_en) = rnd * CPFEM_odd_stress
          CPFEM_dcsde(1:6,1:6,IP,cp_en) = CPFEM_odd_jacobian * math_identity2nd(6)
        else  
          ! translate from P to CS
          Kirchhoff = math_mul33x33(materialpoint_P(1:3,1:3,IP, cp_en), math_transpose3x3(materialpoint_F(1:3,1:3,IP,cp_en)))
          J_inverse  = 1.0_pReal / math_det3x3(materialpoint_F(1:3,1:3,IP,cp_en))
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
    
    ! --+>> COLLECTION OF FEM INPUT WITH RETURNING OF RANDOMIZED ODD STRESS AND JACOBIAN <<+-- 
    case (3,4,5)
      if (mode == 4) then
        CPFEM_dcsde_knownGood = CPFEM_dcsde  ! --+>> BACKUP JACOBIAN FROM FORMER CONVERGED INC
      else if (mode == 5) then
        CPFEM_dcsde = CPFEM_dcsde_knownGood  ! --+>> RESTORE CONSISTENT JACOBIAN FROM FORMER CONVERGED INC
      end if
      call random_number(rnd)
      rnd = 2.0_pReal * rnd - 1.0_pReal
      materialpoint_Temperature(IP,cp_en) = Temperature
      materialpoint_F0(1:3,1:3,IP,cp_en) = ffn
      materialpoint_F(1:3,1:3,IP,cp_en) = ffn1
      CPFEM_cs(1:6,IP,cp_en) = rnd * CPFEM_odd_stress
      CPFEM_dcsde(1:6,1:6,IP,cp_en) = CPFEM_odd_jacobian * math_identity2nd(6)
      CPFEM_calc_done = .false.
    
    ! --+>> RECYCLING OF FORMER RESULTS (MARC SPECIALTY) <<+--
    case (6)
      ! do nothing
    ! --+>> RESTORE CONSISTENT JACOBIAN FROM FORMER CONVERGED INC
    case (7)
      CPFEM_dcsde = CPFEM_dcsde_knownGood
    
  end select

  ! return the local stress and the jacobian from storage
  cauchyStress = CPFEM_cs(1:6,IP,cp_en)
  jacobian = CPFEM_dcsdE(1:6,1:6,IP,cp_en)
  
  ! copy P and dPdF to the output variables 
  pstress = materialpoint_P(1:3,1:3,IP,cp_en)
  dPdF = materialpoint_dPdF(1:3,1:3,1:3,1:3,IP,cp_en)

  ! warning for zero stiffness
  if (all(abs(jacobian) < 1e-10_pReal)) then
    call IO_warning(601,cp_en,IP)
  endif
  
  if (mode < 6 .and. debug_verbosity > 0 .and. ((debug_e == cp_en .and. debug_i == IP) .or. .not. debug_selectiveDebugger)) then
    !$OMP CRITICAL (write2out)
      write(6,'(a,i5,x,i2,/,12(x),6(f10.3,x)/)') '<< CPFEM >> stress/MPa at el ip ', cp_en, IP, cauchyStress/1e6
      write(6,'(a,i5,x,i2,/,6(12(x),6(f10.3,x)/))') '<< CPFEM >> jacobian/GPa at el ip ', cp_en, IP, transpose(jacobian)/1e9
      call flush(6)
    !$OMP END CRITICAL (write2out)
  endif
  
  ! remember extreme values of stress and jacobian
  if (mode < 3) then
    cauchyStress33 = math_Mandel6to33(cauchyStress)
    if (maxval(cauchyStress33) > debug_stressMax) then                        
      debug_stressMaxLocation = (/cp_en, IP/)
      debug_stressMax = maxval(cauchyStress33)
    endif
    if (minval(cauchyStress33) < debug_stressMin) then
      debug_stressMinLocation = (/cp_en, IP/)
      debug_stressMin = minval(cauchyStress33)
    endif
    jacobian3333 = math_Mandel66to3333(jacobian)
    if (maxval(jacobian3333) > debug_jacobianMax) then
      debug_jacobianMaxLocation = (/cp_en, IP/)
      debug_jacobianMax = maxval(jacobian3333)
    endif
    if (minval(jacobian3333) < debug_jacobianMin) then
      debug_jacobianMinLocation = (/cp_en, IP/)
      debug_jacobianMin = minval(jacobian3333)
    endif
  endif
  
  ! return temperature
  if (theTime > 0.0_pReal) Temperature = materialpoint_Temperature(IP,cp_en)  ! homogenized result except for potentially non-isothermal starting condition.

end subroutine

END MODULE CPFEM
