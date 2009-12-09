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

logical ::                                        CPFEM_init_done     = .false., &    ! remember whether init has been done already
                                                  CPFEM_calc_done     = .false.       ! remember whether first IP has already calced the results
 

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
subroutine CPFEM_general(mode, ffn, ffn1, Temperature, dt, element, IP, cauchyStress, jacobian, ngens)
  ! note: cauchyStress = Cauchy stress cs(6) and jacobian = Consistent tangent dcs/dE

  !*** variables and functions from other modules ***!
  use prec, only:                                     pReal, &
                                                      pInt, &
                                                      prec_init
  use numerics, only:                                 numerics_init, & 
                                                      relevantStrain, &
                                                      iJacoStiffness
  use debug, only:                                    debug_init
  use FEsolving, only:                                FE_init, &
                                                      parallelExecution, &
                                                      outdatedFFN1, &
                                                      terminallyIll, &
                                                      cycleCounter, &
                                                      theInc, &
                                                      theTime, &
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
                                                      mesh_maxNips
  use lattice, only:                                  lattice_init
  use material, only:                                 material_init, &
                                                      homogenization_maxNgrains
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
                                                      materialpoint_Temperature, &
                                                      materialpoint_stressAndItsTangent, &
                                                      materialpoint_postResults
  use IO, only:                                       IO_init
  use cpfem_interface
  
  implicit none
  
  !*** input variables ***!
  integer(pInt), intent(in) ::                        element, &          ! FE element number
                                                      IP, &               ! FE integration point number
                                                      ngens               ! size of stress strain law
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
  real(pReal), dimension(ngens), intent(out) ::       cauchyStress        ! stress vector in Mandel notation
  real(pReal), dimension(ngens,ngens), intent(out) :: jacobian            ! jacobian in Mandel notation
  
  !*** local variables ***!
  real(pReal)                                         J_inverse, &        ! inverse of Jacobian
                                                      rnd
  real(pReal), dimension (3,3) ::                     Kirchhoff
  real(pReal), dimension (3,3,3,3) ::                 H, &
                                                      H_sym
  integer(pInt)                                       cp_en, &            ! crystal plasticity element number
                                                      i, &
                                                      j, &
                                                      k, &
                                                      l, &
                                                      m, &
                                                      n
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
    call mpie_cpfem_init()
    CPFEM_init_done = .true.
  endif

  cp_en = mesh_FEasCP('elem',element)
  
  if (cp_en == 1 .and. IP == 1) then
  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,*) '#####################################'
    write(6,'(a10,x,f8.4,x,a10,x,i4,x,a10,x,i3,x,a16,x,i2,x,a16,x,i2)') &
    'theTime',theTime,'theInc',theInc,'cycleCounter',cycleCounter,'computationMode',mode
    write(6,*) '#####################################'
    call flush (6)
  !$OMP END CRITICAL (write2out)
  endif

  ! according to our "mode" we decide what to do
  select case (mode)
    
    ! --+>> REGULAR COMPUTATION (WITH AGING OF RESULTS IF MODE == 1) <<+-- 
    case (1,2)
      ! age results if mode == 1
      if (mode == 1) then
        crystallite_F0  = crystallite_partionedF                          ! crystallite deformation (_subF is perturbed...)
        crystallite_Fp0 = crystallite_Fp                                  ! crystallite plastic deformation
        crystallite_Lp0 = crystallite_Lp                                  ! crystallite plastic velocity
        crystallite_Tstar0_v = crystallite_Tstar_v                        ! crystallite 2nd Piola Kirchhoff stress 
        forall ( i = 1:homogenization_maxNgrains, &
                 j = 1:mesh_maxNips, &
                 k = 1:mesh_NcpElems ) &
          constitutive_state0(i,j,k)%p = constitutive_state(i,j,k)%p      ! microstructure of crystallites
  !$OMP CRITICAL (write2out)
        write(6,'(a10,/,4(3(e20.8,x),/))') 'aged state',constitutive_state(1,1,1)%p
  !$OMP END CRITICAL (write2out)
        do k = 1,mesh_NcpElems
          do j = 1,mesh_maxNips
            if (homogenization_sizeState(j,k) > 0_pInt) &
              homogenization_state0(j,k)%p = homogenization_state(j,k)%p  ! internal state of homogenization scheme
          enddo
        enddo
      endif

      ! deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
      if (terminallyIll .or. outdatedFFN1 .or. any(abs(ffn1 - materialpoint_F(:,:,IP,cp_en)) > relevantStrain)) then
        if (.not. terminallyIll .and. .not. outdatedFFN1) then 
  !$OMP CRITICAL (write2out)
          write(6,'(a11,x,i5,x,i2,x,a10,/,3(3(f10.6,x),/))') 'outdated at',cp_en,IP,'FFN1 now:',ffn1(:,1),ffn1(:,2),ffn1(:,3)
  !$OMP END CRITICAL (write2out)
          outdatedFFN1 = .true.
        endif
        CPFEM_cs(1:ngens,IP,cp_en)              = CPFEM_odd_stress
        CPFEM_dcsde(1:ngens,1:ngens,IP,cp_en)   = CPFEM_odd_jacobian*math_identity2nd(ngens)
      
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
          CPFEM_calc_done = .true.
        endif
        
        if (terminallyIll) then
          CPFEM_cs(1:ngens,IP,cp_en)              = CPFEM_odd_stress
          CPFEM_dcsde(1:ngens,1:ngens,IP,cp_en)   = CPFEM_odd_jacobian*math_identity2nd(ngens)
        else  
        !  translate from P to CS
          Kirchhoff = math_mul33x33(materialpoint_P(:,:,IP, cp_en),transpose(materialpoint_F(:,:,IP, cp_en)))
          J_inverse  = 1.0_pReal/math_det3x3(materialpoint_F(:,:,IP, cp_en))
          CPFEM_cs(1:ngens,IP,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff)

        !  translate from dP/dF to dCS/dE
          H = 0.0_pReal
          forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
            H(i,j,k,l) = H(i,j,k,l) + &
                          materialpoint_F(j,m,IP,cp_en) * &
                          materialpoint_F(l,n,IP,cp_en) * &
                          materialpoint_dPdF(i,m,k,n,IP,cp_en) - &
                          math_I3(j,l)*materialpoint_F(i,m,IP,cp_en)*materialpoint_P(k,m,IP,cp_en) + &
                          0.5_pReal*(math_I3(i,k)*Kirchhoff(j,l) + math_I3(j,l)*Kirchhoff(i,k) + &
                                     math_I3(i,l)*Kirchhoff(j,k) + math_I3(j,k)*Kirchhoff(i,l))
          forall(i=1:3,j=1:3,k=1:3,l=1:3) &
            H_sym(i,j,k,l)= 0.25_pReal*(H(i,j,k,l)+H(j,i,k,l)+H(i,j,l,k)+H(j,i,l,k))  ! where to use the symmetric version??
          CPFEM_dcsde(1:ngens,1:ngens,IP,cp_en) = math_Mandel3333to66(J_inverse*H)
        endif
      endif
    
    ! --+>> COLLECTION OF FEM INPUT WITH RETURNING OF ODD STRESS AND JACOBIAN <<+-- 
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
      CPFEM_cs(1:ngens,IP,cp_en)            = rnd*CPFEM_odd_stress
      CPFEM_dcsde(1:ngens,1:ngens,IP,cp_en) = CPFEM_odd_jacobian*math_identity2nd(ngens)
      CPFEM_calc_done = .false.
    
    ! --+>> RECYCLING OF FORMER RESULTS (MARC SPECIALTY) <<+--
    case (6,7)
      if (mode == 7) CPFEM_dcsde = CPFEM_dcsde_knownGood  ! --+>> RESTORE CONSISTENT JACOBIAN FROM FORMER CONVERGED INC
      ! do nothing
    
  end select

  ! return the local stress and the jacobian from storage
  cauchyStress(1:ngens) = CPFEM_cs(1:ngens,IP,cp_en)
  jacobian(1:ngens,1:ngens) = CPFEM_dcsdE(1:ngens,1:ngens,IP,cp_en)
  if (IP == 1 .and. cp_en == 1) then
  !$OMP CRITICAL (write2out)
    write(6,'(a,/,6(6(f10.3,x)/))') 'jacobian/GPa at ip 1 el 1',jacobian/1e9
    call flush(6)
  !$OMP END CRITICAL (write2out)
  endif
  
  ! return temperature
  if (theInc > 0_pInt) Temperature = materialpoint_Temperature(IP,cp_en)  ! homogenized result except for potentially non-isothermal starting condition.
  return

end subroutine

END MODULE CPFEM
 
