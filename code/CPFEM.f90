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
  allocate(CPFEM_dcsdE(6,6,mesh_maxNips,mesh_NcpElems)) ;           CPFEM_dcsde           = 0.0_pReal
  allocate(CPFEM_dcsdE_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dcsde_knownGood = 0.0_pReal

  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,*) '<<<+-  cpfem init  -+>>>'
    write(6,*)
    write(6,'(a32,x,6(i5,x))') 'CPFEM_cs:              ', shape(CPFEM_cs)
    write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsde:           ', shape(CPFEM_dcsde)
    write(6,'(a32,x,6(i5,x))') 'CPFEM_dcsde_knownGood: ', shape(CPFEM_dcsde_knownGood)
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
  ! note: cauchyStress = Cauchy stress cs(6) and jacobian = Consistent tangent dcs/de

  !*** variables and functions from other modules ***!
  use prec, only:                                     pReal, &
                                                      pInt
  use numerics, only:                                 numerics_init, & 
                                                      relevantStrain, &
                                                      iJacoStiffness
  use debug, only:                                    debug_init
  use FEsolving, only:                                FE_init, &
                                                      parallelExecution, &
                                                      outdatedFFN1, &
                                                      cycleCounter, &
                                                      theInc, &
                                                      theCycle, &
                                                      theLovl, &
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
                                                      mesh_maxNips, &
													  mesh_element
  use lattice, only:                                  lattice_init
  use material, only:                                 material_init, &
                                                      homogenization_maxNgrains, &
													  homogenization_Ngrains
  use constitutive, only:                             constitutive_init,&
                                                      constitutive_state0,constitutive_state
  use crystallite, only:                              crystallite_init, &
                                                      crystallite_F0, &
                                                      crystallite_partionedF, &
                                                      crystallite_Fp0, &
                                                      crystallite_Fp, &
                                                      crystallite_Lp0, &
                                                      crystallite_Lp
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
  
  implicit none
  
  !*** input variables ***!
  integer(pInt), intent(in) ::                        element, &          ! FE element number
                                                      IP, &               ! FE integration point number
                                                      ngens               ! size of stress strain law
  real(pReal), intent(inout) ::                       Temperature         ! temperature
  real(pReal), intent(in) ::                          dt                  ! time increment
  real(pReal), dimension (3,3), intent(in) ::         ffn, &              ! deformation gradient for t=t0
                                                      ffn1                ! deformation gradient for t=t1
  integer(pInt), intent(in) ::                        mode                ! computation mode  1: regular computation with aged results
                                                                          !                   2: regular computation
                                                                          !                   3: collection of FEM data
                                                                          !                   4: recycling of former results (MARC speciality)
                                                                          !                   5: record tangent from former converged inc
                                                                          !                   6: restore tangent from former converged inc
  
  !*** output variables ***!
  real(pReal), dimension(ngens), intent(out) ::       cauchyStress        ! stress vector in Mandel notation
  real(pReal), dimension(ngens,ngens), intent(out) :: jacobian            ! jacobian in Mandel notation
  
  !*** local variables ***!
  real(pReal)                                         J_inverse           ! inverse of Jacobian
  real(pReal), dimension (3,3) ::                     Kirchhoff
  real(pReal), dimension (3,3,3,3) ::                 H, &
                                                      H_sym
  integer(pInt)                                       cp_en, &            ! crystal plasticity element number
                                                      e, &
                                                      g, &
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
    call numerics_init()
    call debug_init()
    call math_init()
    call FE_init()
    call mesh_init()
    call lattice_init()
    call material_init()
    call constitutive_init()
    call crystallite_init(Temperature)         ! (have to) use temperature of first IP for whole model
    call homogenization_init(Temperature)
    call CPFEM_init()
    CPFEM_init_done = .true.
  endif

  cp_en = mesh_FEasCP('elem',element)
  
  if (cp_en == 1 .and. IP == 1) then
    write(6,*) '#####################################'
    write(6,'(a10,1x,f8.4,1x,a10,1x,i4,1x,a10,1x,i3,1x,a10,1x,i2,x,a10,1x,i2)') &
    'theTime',theTime,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,&
    'mode',mode
    write(6,*) '#####################################'
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
        do e = 1,mesh_NcpElems		
		  do g = 1,homogenization_Ngrains(mesh_element(3,e))
		    do i = 1,mesh_maxNips
               constitutive_state0(g,i,e)%p = constitutive_state(g,i,e)%p      ! microstructure of crystallites
            enddo
		  enddo
		enddo      
		write(6,'(a10,/,4(3(f10.3,x),/))') 'aged state',constitutive_state(1,1,1)%p/1e6
        do e = 1,mesh_NcpElems
          do i = 1,mesh_maxNips
            if (homogenization_sizeState(i,e) > 0_pInt) &
              homogenization_state0(i,e)%p = homogenization_state(i,e)%p  ! internal state of homogenization scheme
          enddo
        enddo
      endif

      ! deformation gradient outdated or any actual deformation gradient differs more than relevantStrain from the stored one
      if (outdatedFFN1 .or. any(abs(ffn1 - materialpoint_F(:,:,IP,cp_en)) > relevantStrain)) then
        if (.not. outdatedFFN1) & 
          write(6,'(a11,x,i5,x,i2,x,a10,/,3(3(f10.3,x),/))') 'outdated at',cp_en,IP,'FFN1 now:',ffn1(:,1),ffn1(:,2),ffn1(:,3)
        outdatedFFN1 = .true.
        CPFEM_cs(1:ngens,IP,cp_en)              = CPFEM_odd_stress
        CPFEM_dcsde(1:ngens,1:ngens,IP,cp_en)   = CPFEM_odd_jacobian*math_identity2nd(ngens)
      
      ! deformation gradient is not outdated
      else
        ! set flag for Jacobian update
        updateJaco = (mod(cycleCounter-4,4_pInt*iJacoStiffness)==0)
        
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
    
    ! --+>> COLLECTION OF FEM DATA AND RETURN OF ODD STRESS AND JACOBIAN <<+-- 
    case (3)
	  if (IP==1.AND.cp_en==1) write(6,*) 'Temp from CPFEM', Temperature
      materialpoint_Temperature(IP,cp_en)   = Temperature
      materialpoint_F0(:,:,IP,cp_en)        = ffn
      materialpoint_F(:,:,IP,cp_en)         = ffn1
      CPFEM_cs(1:ngens,IP,cp_en)            = CPFEM_odd_stress
      CPFEM_dcsde(1:ngens,1:ngens,IP,cp_en) = CPFEM_odd_jacobian*math_identity2nd(ngens)
      CPFEM_calc_done = .false.
    
    ! --+>> RECYCLING OF FORMER RESULTS (MARC SPECIALTY) <<+--
    case (4)
      ! do nothing
    
    ! --+>> RECORD JACOBIAN FROM FORMER CONVERGED INC <<+--
    case (5)
       CPFEM_dcsde_knownGood = CPFEM_dcsde
    
    ! --+>> RESTORE CONSISTENT JACOBIAN FROM FORMER CONVERGED INC <<+--
    case (6)
       CPFEM_dcsde = CPFEM_dcsde_knownGood
       
  end select

  ! return the local stress and the jacobian from storage
  cauchyStress(1:ngens) = CPFEM_cs(1:ngens,IP,cp_en)
  jacobian(1:ngens,1:ngens) = CPFEM_dcsdE(1:ngens,1:ngens,IP,cp_en)
  ! return temperature
  if (theInc > 0_pInt) Temperature = materialpoint_Temperature(IP,cp_en)  ! homogenized result except for potentially non-isothermal starting condition.
  
  return

end subroutine



 END MODULE CPFEM
 