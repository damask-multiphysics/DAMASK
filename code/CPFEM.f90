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
!***    allocate the arrays defined in module CPFEM    ***
!***    and initialize them                            ***
!*********************************************************
subroutine CPFEM_init()

  use prec, only:                                 pInt
  use FEsolving, only:                            parallelExecution, &
                                                  symmetricSolver
  use mesh, only:                                 mesh_NcpElems, &
                                                  mesh_maxNips

  implicit none

  ! initialize stress and jacobian to zero 
  allocate(CPFEM_cs(6,mesh_maxNips,mesh_NcpElems)) ;                CPFEM_cs              = 0.0_pReal
  allocate(CPFEM_dcsdE(6,6,mesh_maxNips,mesh_NcpElems)) ;           CPFEM_dcsdE           = 0.0_pReal
  allocate(CPFEM_dcsdE_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dcsdE_knownGood = 0.0_pReal

  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,*) '<<<+-  cpfem init  -+>>>'
    write(6,*) '$Id$'
    write(6,*)
    write(6,'(a32,x,6(i5,x))') 'CPFEM_cs:              ', shape(CPFEM_cs)
    write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
    write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
    write(6,*)
    write(6,*) 'parallelExecution:    ', parallelExecution
    write(6,*) 'symmetricSolver:      ', symmetricSolver
    call flush(6)
  !$OMP END CRITICAL (write2out)
  return

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
                                                      pInt, &
                                                      prec_init
  use numerics, only:                                 numerics_init, & 
                                                      relevantStrain, &
                                                      defgradTolerance, &
                                                      iJacoStiffness
  use debug, only:                                    debug_init, &
                                                      debug_g, &
                                                      debug_i, &
                                                      debug_e, &
                                                      debugger, &
                                                      selectiveDebugger, &
                                                      verboseDebugger
  use FEsolving, only:                                FE_init, &
                                                      parallelExecution, &
                                                      outdatedFFN1, &
                                                      terminallyIll, &
                                                      cycleCounter, &
                                                      theInc, &
                                                      theTime, &
                                                      theDelta, &
                                                      FEsolving_execElem, &
                                                      FEsolving_execIP
  use math, only:                                     math_init, &
                                                      math_identity2nd, &
                                                      math_mul33x33, &
                                                      math_det3x3, &
                                                      math_I3, &
                                                      math_Mandel3333to66, &
                                                      math_Mandel33to6
  use mesh, only:                                     mesh_init, &
                                                      mesh_FEasCP, &
                                                      mesh_NcpElems, &
                                                      mesh_maxNips, &
                                                      mesh_element, &
                                                      FE_Nips
  use lattice, only:                                  lattice_init
  use material, only:                                 material_init, &
                                                      homogenization_maxNgrains, &
                                                      microstructure_elemhomo
  use constitutive, only:                             constitutive_init,&
                                                      constitutive_state0,constitutive_state
  use crystallite, only:                              crystallite_init, &
                                                      crystallite_F0, &
                                                      crystallite_partionedF, &
                                                      crystallite_Fp0, &
                                                      crystallite_Fp, &
                                                      crystallite_Lp0, &
                                                      crystallite_Lp, &
                                                      crystallite_Tstar0_v, &
                                                      crystallite_Tstar_v
  use homogenization, only:                           homogenization_init, &
                                                      homogenization_sizeState, &
                                                      homogenization_state, &
                                                      homogenization_state0, &
                                                      materialpoint_F, &
                                                      materialpoint_F0, &
                                                      materialpoint_P, &
                                                      materialpoint_dPdF, &
                                                      materialpoint_results, &
                                                      materialpoint_Temperature, &
                                                      materialpoint_stressAndItsTangent, &
                                                      materialpoint_postResults
  use IO, only:                                       IO_init
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
  real(pReal), dimension (3,3), intent(out) ::        pstress                   ! Piola-Kirchhoff stress in Matrix notation
  real(pReal), dimension (3,3,3,3), intent(out) ::    dPdF                ! 
                 
  !*** local variables ***!
  real(pReal)                                         J_inverse, &        ! inverse of Jacobian
                                                      rnd
  real(pReal), dimension (3,3) ::                     Kirchhoff           ! Piola-Kirchhoff stress in Matrix notation
  real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                      H
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
  
  cp_en = mesh_FEasCP('elem',element)
  selectiveDebugger = (cp_en == debug_e .and. IP == debug_i)
  
  if (selectiveDebugger) then
    !$OMP CRITICAL (write2out)
      write(6,*)
      write(6,'(a)') '#######################################################'
      write(6,'(a32,x,i5,x,i2)') 'reporting for element, ip:',cp_en,IP
      write(6,'(a32,x,f15.7)') 'theTime',theTime
      write(6,'(a32,x,f15.7)') 'theDelta',theDelta
      write(6,'(a32,x,i8)') 'theInc',theInc
      write(6,'(a32,x,i8)') 'cycleCounter',cycleCounter
      write(6,'(a32,x,i8)') 'computationMode',mode
      write(6,'(a)') '#######################################################'
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
        crystallite_Tstar0_v = crystallite_Tstar_v                        ! crystallite 2nd Piola Kirchhoff stress 
        forall ( i = 1:homogenization_maxNgrains, &
                 j = 1:mesh_maxNips, &
                 k = 1:mesh_NcpElems ) &
          constitutive_state0(i,j,k)%p = constitutive_state(i,j,k)%p      ! microstructure of crystallites
        if (selectiveDebugger) then
          !$OMP CRITICAL (write2out)
            write(6,'(a32,x,i8,x,i2,/,4(3(e20.8,x),/))') '°°° AGED state of grain 1, element ip',cp_en,IP, &
                                                         constitutive_state(1,IP,cp_en)%p
          !$OMP END CRITICAL (write2out)
        endif
        do k = 1,mesh_NcpElems
          do j = 1,mesh_maxNips
            if (homogenization_sizeState(j,k) > 0_pInt) &
              homogenization_state0(j,k)%p = homogenization_state(j,k)%p  ! internal state of homogenization scheme
          enddo
        enddo
      endif

      if (mode == 8 .or. mode == 9) then                                  ! Abaqus explicit skips collect
        materialpoint_Temperature(IP,cp_en)   = Temperature
        materialpoint_F0(:,:,IP,cp_en)        = ffn
        materialpoint_F(:,:,IP,cp_en)         = ffn1
      endif

      ! deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
      if (terminallyIll .or. outdatedFFN1 .or. any(abs(ffn1 - materialpoint_F(:,:,IP,cp_en)) > defgradTolerance)) then
        if (.not. terminallyIll .and. .not. outdatedFFN1) then 
          !$OMP CRITICAL (write2out)
            write(6,'(a32,x,i5,x,i2)') '°°° OUTDATED at element ip',cp_en,IP
            write(6,'(a32,/,3(3(f10.6,x),/))') '    FFN1 now:',ffn1(:,1),ffn1(:,2),ffn1(:,3)
          !$OMP END CRITICAL (write2out)
          outdatedFFN1 = .true.
        endif
        call random_number(rnd)
        if (rnd < 0.5_pReal) rnd = 1.0_pReal - rnd
        CPFEM_cs(:,IP,cp_en)              = rnd*CPFEM_odd_stress
        CPFEM_dcsde(:,:,IP,cp_en)         = CPFEM_odd_jacobian*math_identity2nd(6)
      
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
          call materialpoint_stressAndItsTangent(updateJaco, dt)          ! calculate stress and its tangent
          call materialpoint_postResults(dt)                              ! post results
          
        ! parallel computation and calulation not yet done
        elseif (.not. CPFEM_calc_done) then
          call materialpoint_stressAndItsTangent(updateJaco, dt)          ! calculate stress and its tangent (parallel execution inside)
          call materialpoint_postResults(dt)                              ! post results
          do e = FEsolving_execElem(1),FEsolving_execElem(2)              ! loop over all parallely processed elements
            if (microstructure_elemhomo(mesh_element(4,e))) then          ! dealing with homogeneous element?
              forall (i = 2:FE_Nips(mesh_element(2,e)))                   ! copy results of first IP to all others
                materialpoint_P(:,:,i,e)        = materialpoint_P(:,:,1,e) 
                materialpoint_F(:,:,i,e)        = materialpoint_F(:,:,1,e) 
                materialpoint_dPdF(:,:,:,:,i,e) = materialpoint_dPdF(:,:,:,:,1,e)
                materialpoint_results(:,i,e)    = materialpoint_results(:,1,e)
              end forall
            endif
          enddo
          CPFEM_calc_done = .true.
        endif
        
        if ( terminallyIll ) then
          call random_number(rnd)
          if (rnd < 0.5_pReal) rnd = 1.0_pReal - rnd
          CPFEM_cs(:,IP,cp_en)              = rnd*CPFEM_odd_stress
          CPFEM_dcsde(:,:,IP,cp_en)         = CPFEM_odd_jacobian*math_identity2nd(6)
        else  
        !  translate from P to CS
          Kirchhoff = math_mul33x33(materialpoint_P(:,:,IP, cp_en),transpose(materialpoint_F(:,:,IP, cp_en)))
          J_inverse  = 1.0_pReal/math_det3x3(materialpoint_F(:,:,IP, cp_en))
          CPFEM_cs(:,IP,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff)

        !  translate from dP/dF to dCS/dE
          H = 0.0_pReal
          do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
!          forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
            H(i,j,k,l) = H(i,j,k,l) + &
                          materialpoint_F(j,m,IP,cp_en) * &
                          materialpoint_F(l,n,IP,cp_en) * &
                          materialpoint_dPdF(i,m,k,n,IP,cp_en) - &
                          math_I3(j,l)*materialpoint_F(i,m,IP,cp_en)*materialpoint_P(k,m,IP,cp_en) + &
                          0.5_pReal*(math_I3(i,k)*Kirchhoff(j,l) + math_I3(j,l)*Kirchhoff(i,k) + &
                                     math_I3(i,l)*Kirchhoff(j,k) + math_I3(j,k)*Kirchhoff(i,l))
          enddo; enddo; enddo; enddo; enddo; enddo
!          forall(i=1:3,j=1:3,k=1:3,l=1:3) &
!            H_sym(i,j,k,l) = 0.25_pReal*(H(i,j,k,l)+H(j,i,k,l)+H(i,j,l,k)+H(j,i,l,k))  ! where to use this symmetric version ??
          CPFEM_dcsde(:,:,IP,cp_en) = math_Mandel3333to66(J_inverse*H)                 ! should this use the symmetrized H ??
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
      if (rnd < 0.5_pReal) rnd = 1.0_pReal - rnd
      materialpoint_Temperature(IP,cp_en)   = Temperature
      materialpoint_F0(:,:,IP,cp_en)        = ffn
      materialpoint_F(:,:,IP,cp_en)         = ffn1
      CPFEM_cs(:,IP,cp_en)                  = rnd*CPFEM_odd_stress
      CPFEM_dcsde(:,:,IP,cp_en)             = CPFEM_odd_jacobian*math_identity2nd(6)
      CPFEM_calc_done = .false.
    
    ! --+>> RECYCLING OF FORMER RESULTS (MARC SPECIALTY) <<+--
    case (6)
      ! do nothing
    ! --+>> RESTORE CONSISTENT JACOBIAN FROM FORMER CONVERGED INC
    case (7)
      CPFEM_dcsde = CPFEM_dcsde_knownGood
    
  end select

  ! return the local stress and the jacobian from storage
  cauchyStress(:) = CPFEM_cs(:,IP,cp_en)
  jacobian(:,:)   = CPFEM_dcsdE(:,:,IP,cp_en)
  
  ! copy P and dPdF to the output variables 
  pstress(:,:)   = materialpoint_P(:,:,IP,cp_en)
  dPdF(:,:,:,:)  = materialpoint_dPdF(:,:,:,:,IP,cp_en)
  if ((debugger .and. selectiveDebugger) .and. &
      mode < 6) then
    !$OMP CRITICAL (write2out)
      write(6,'(a,x,i2,x,a,x,i4,/,6(f10.3,x)/)') 'stress/MPa at ip', IP, 'el', cp_en, cauchyStress/1e6
      write(6,'(a,x,i2,x,a,x,i4,/,6(6(f10.3,x)/))') 'jacobian/GPa at ip', IP, 'el', cp_en, jacobian/1e9
      call flush(6)
    !$OMP END CRITICAL (write2out)
  endif
  
  ! return temperature
  if (theTime > 0.0_pReal) Temperature = materialpoint_Temperature(IP,cp_en)  ! homogenized result except for potentially non-isothermal starting condition.
  return

end subroutine

END MODULE CPFEM
