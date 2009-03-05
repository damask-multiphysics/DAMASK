!##############################################################
 MODULE CPFEM
!##############################################################
!    *** CPFEM engine ***
!
 use prec, only: pReal,pInt
 implicit none
!
! ****************************************************************
! *** General variables for the material behaviour calculation ***
! ****************************************************************
 real(pReal), dimension (:,:),          allocatable :: CPFEM_Temperature
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_ffn_bar          !average FFN per IP
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_ffn              !individual FFN per grain
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_ffn1_bar         !average FFN1 per IP
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_ffn1             !individual FFN1 per grain
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_PK1_bar          !average PK1 per IP
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_PK1              !individual PK1 per grain
 real(pReal), dimension (:,:,:,:,:,:),  allocatable :: CPFEM_dPdF_bar         !average dPdF per IP
 real(pReal), dimension (:,:,:,:,:,:),  allocatable :: CPFEM_dPdF_bar_old     !old average dPdF per IP
 real(pReal), dimension (:,:,:,:,:,:,:),allocatable :: CPFEM_dPdF             !individual dPdF per grain
 real(pReal), dimension (:,:,:),        allocatable :: CPFEM_stress_bar
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_jaco_bar
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_jaco_knownGood
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_results
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_Lp_old
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_Lp_new
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_Fp_old
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_Fp_new
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_Fe_new
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_Tstar_v

 logical, dimension (:,:,:), allocatable :: crystallite_converged            !individual convergence flag per grain

 integer(pInt), dimension(:,:), allocatable :: CPFEM_execution_IP
 integer(pInt), dimension(2) :: CPFEM_execution_elem

 integer(pInt) :: CPFEM_Nresults   = 5_pInt    ! phase, volfrac, three Euler angles
 logical :: CPFEM_init_done    = .false.       ! remember whether init has been done already
 logical :: CPFEM_calc_done    = .false.       ! remember whether first IP has already calced the results

 real(pReal), parameter :: CPFEM_odd_stress = 1e15_pReal, CPFEM_odd_jacobian = 1e50_pReal
!
 CONTAINS
!
!*********************************************************
!***    allocate the arrays defined in module CPFEM    ***
!***    and initialize them                            ***
!*********************************************************
 SUBROUTINE CPFEM_init(Temperature)
!
 use prec
 use math, only: math_EulertoR, math_I3, math_identity2nd
 use FEsolving, only: parallelExecution
 use mesh
 use material
 use constitutive
!
 implicit none
!
 real(pReal) Temperature
 integer(pInt) e,i,g
!
!    *** mpie.marc parameters ***
 allocate(CPFEM_Temperature(mesh_maxNips,mesh_NcpElems)) ; CPFEM_Temperature = Temperature
 allocate(CPFEM_ffn_bar(3,3,mesh_maxNips,mesh_NcpElems))
 forall(e=1:mesh_NcpElems,i=1:mesh_maxNips) CPFEM_ffn_bar(:,:,i,e)             = math_I3
 allocate(CPFEM_ffn(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall(g=1:homogenization_maxNgrains,e=1:mesh_NcpElems,i=1:mesh_maxNips) CPFEM_ffn(:,:,g,i,e) = math_I3
 allocate(CPFEM_ffn1_bar(3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1_bar   = CPFEM_ffn_bar
 allocate(CPFEM_ffn1(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1 = CPFEM_ffn
 allocate(CPFEM_PK1_bar(3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_PK1_bar    = 0.0_pReal
 allocate(CPFEM_PK1(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_PK1 = 0.0_pReal
 allocate(CPFEM_dPdF_bar(3,3,3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF_bar = 0.0_pReal
 allocate(CPFEM_dPdF_bar_old(3,3,3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF_bar_old = 0.0_pReal
 allocate(CPFEM_dPdF(3,3,3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF = 0.0_pReal
 allocate(CPFEM_stress_bar(6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_stress_bar = 0.0_pReal
 allocate(CPFEM_jaco_bar(6,6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_jaco_bar   = 0.0_pReal
 allocate(CPFEM_jaco_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_jaco_knownGood = 0.0_pReal
!
!    *** User defined results ***
 allocate(CPFEM_results(CPFEM_Nresults+constitutive_maxSizePostResults,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_results = 0.0_pReal
!
!    *** Plastic velocity gradient ***
 allocate(CPFEM_Lp_old(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Lp_old = 0.0_pReal
 allocate(CPFEM_Lp_new(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Lp_new = 0.0_pReal

!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***
 allocate(CPFEM_Fp_new(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fp_new = 0.0_pReal
 allocate(CPFEM_Fp_old(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall (e=1:mesh_NcpElems,i=1:mesh_maxNips,g=1:homogenization_maxNgrains) &
   CPFEM_Fp_old(:,:,g,i,e) = math_EulerToR(material_EulerAngles(:,g,i,e))  ! plastic def gradient reflects init orientation
!    *** Elastic deformation gradient at (t=t1) ***
 allocate(CPFEM_Fe_new(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fe_new = 0.0_pReal
!    *** Stress vector at (t=t1) ***
 allocate(CPFEM_Tstar_v(6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Tstar_v = 0.0_pReal
!
 allocate(crystallite_converged(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)); crystallite_converged = .false.

 allocate(CPFEM_execution_IP(2,mesh_NcpElems)); CPFEM_execution_IP = 1_pInt
 forall (e = 1:mesh_NcpElems) CPFEM_execution_IP(2,e) = FE_Nips(mesh_element(2,e))
 CPFEM_execution_elem = (/1,mesh_NcpElems/)

!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) 'CPFEM Initialization'
 write(6,*)
 write(6,*) 'CPFEM_Temperature:    ', shape(CPFEM_Temperature)
 write(6,*) 'CPFEM_ffn_bar:        ', shape(CPFEM_ffn_bar)
 write(6,*) 'CPFEM_ffn:            ', shape(CPFEM_ffn)
 write(6,*) 'CPFEM_ffn1_bar:       ', shape(CPFEM_ffn1_bar)
 write(6,*) 'CPFEM_ffn1:           ', shape(CPFEM_ffn1)
 write(6,*) 'CPFEM_PK1_bar:        ', shape(CPFEM_PK1_bar)
 write(6,*) 'CPFEM_PK1:            ', shape(CPFEM_PK1)
 write(6,*) 'CPFEM_dPdF_bar:       ', shape(CPFEM_dPdF_bar)
 write(6,*) 'CPFEM_dPdF_bar_old:   ', shape(CPFEM_dPdF_bar_old)
 write(6,*) 'CPFEM_dPdF:           ', shape(CPFEM_dPdF)
 write(6,*) 'CPFEM_stress_bar:     ', shape(CPFEM_stress_bar)
 write(6,*) 'CPFEM_jaco_bar:       ', shape(CPFEM_jaco_bar)
 write(6,*) 'CPFEM_jaco_knownGood: ', shape(CPFEM_jaco_knownGood)
 write(6,*) 'CPFEM_results:        ', shape(CPFEM_results)
 write(6,*) 'CPFEM_Lp_old:         ', shape(CPFEM_Lp_old)
 write(6,*) 'CPFEM_Lp_new:         ', shape(CPFEM_Lp_new)
 write(6,*) 'CPFEM_Fp_old:         ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:         ', shape(CPFEM_Fp_new)
 write(6,*) 'CPFEM_Fe_new:         ', shape(CPFEM_Fe_new)
 write(6,*) 'CPFEM_Tstar_v:        ', shape(CPFEM_Tstar_v)
 write(6,*) 'crystallite_converged:', shape(crystallite_converged)
 write(6,*)
 write(6,*) 'parallelExecution:    ', parallelExecution
 call flush(6)
!$OMP END CRITICAL (write2out)
 return
!
 END SUBROUTINE
!
!
!***********************************************************************
!***    perform initialization at first call, update variables and   ***
!***    call the actual material model                               ***
!
!     CPFEM_mode             computation mode (regular, collection, recycle)
!     ffn                    deformation gradient for t=t0
!     ffn1                   deformation gradient for t=t1
!     Temperature            temperature
!     CPFEM_dt               time increment
!     CPFEM_en               element number
!     CPFEM_in               intergration point number
!     CPFEM_stress           stress vector in Mandel notation
!     CPFEM_updateJaco       flag to initiate computation of Jacobian
!     CPFEM_jaco             jacobian in Mandel notation
!     CPFEM_ngens            size of stress strain law
!***********************************************************************
 SUBROUTINE CPFEM_general(CPFEM_mode, ffn, ffn1, Temperature, CPFEM_dt,&
                          CPFEM_en, CPFEM_in, CPFEM_stress, CPFEM_updateJaco, CPFEM_jaco, CPFEM_ngens)
! note: CPFEM_stress = Cauchy stress cs(6) and CPFEM_jaco = Consistent tangent dcs/de
!
 use prec, only: pReal,pInt
 use FEsolving
 use debug
 use math
 use mesh, only: mesh_init,mesh_FEasCP, mesh_NcpElems, mesh_maxNips, mesh_element
 use lattice, only: lattice_init
 use material
 use constitutive, only: constitutive_init,constitutive_state_old,constitutive_state_new
 implicit none
!
 integer(pInt) CPFEM_en, CPFEM_in, cp_en, CPFEM_ngens, i,j,k,l,m,n
 real(pReal), dimension (3,3)        :: ffn,ffn1,Kirchhoff_bar
 real(pReal), dimension (3,3,3,3)    :: H_bar, H_bar_sym
 real(pReal), dimension(CPFEM_ngens) :: CPFEM_stress
 real(pReal), dimension(CPFEM_ngens,CPFEM_ngens) :: CPFEM_jaco
 real(pReal) Temperature,CPFEM_dt,J_inverse
 integer(pInt) CPFEM_mode               ! 1: regular computation with aged results&
                                        ! 2: regular computation&
                                        ! 3: collection of FEM data&
                                        ! 4: recycling of former results (MARC speciality)&
                                        ! 5: record tangent from former converged inc&
                                        ! 6: restore tangent from former converged inc
 logical       CPFEM_updateJaco
!
 if (.not. CPFEM_init_done) then        ! initialization step (three dimensional stress state check missing?)
   call math_init()
   call FE_init()
   call mesh_init()
   call lattice_init()
   call material_init()
   call constitutive_init()
   write (6,*) 'call CPFEM init'
   call CPFEM_init(Temperature)
   CPFEM_init_done = .true.
 endif
!
 cp_en = mesh_FEasCP('elem',CPFEM_en)
 if (cp_en == 1 .and. CPFEM_in == 1) then
    write(6,'(a10,1x,f8.4,1x,a10,1x,i4,1x,a10,1x,i3,1x,a10,1x,i2,x,a10,1x,i2)') &
    'theTime',theTime,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,&
    'mode',CPFEM_mode
 endif
!
 select case (CPFEM_mode)
    case (1,2)     ! regular computation (with aging of results if mode == 1)
       if (CPFEM_mode == 1) then                  ! age results at start of new increment
         CPFEM_Lp_old            = CPFEM_Lp_new
         CPFEM_Fp_old            = CPFEM_Fp_new
         forall (i = 1:homogenization_maxNgrains,&
                 j = 1:mesh_maxNips, &
                 k = 1:mesh_NcpElems) &
           constitutive_state_old(i,j,k)%p = constitutive_state_new(i,j,k)%p
         write (6,*) 'results aged.'
       endif

       if (outdatedFFN1 .or. any(abs(ffn1 - CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en)) > relevantStrain)) then
         if (.not. outdatedFFN1) write(6,'(i5,x,i2,x,a10,/,3(3(f10.3,x),/))') cp_en,CPFEM_in,'FFN1 now:',ffn1(:,1),ffn1(:,2),ffn1(:,3)
         outdatedFFN1 = .true.
         CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_stress
         CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       else
         if (.not. parallelExecution) then
           CPFEM_execution_elem(1) = cp_en
           CPFEM_execution_elem(2) = cp_en
           CPFEM_execution_IP(1,cp_en) = CPFEM_in
           CPFEM_execution_IP(2,cp_en) = CPFEM_in
           call CPFEM_MaterialPoint(CPFEM_updateJaco, CPFEM_dt)
         elseif (.not. CPFEM_calc_done) then
           call CPFEM_MaterialPoint(CPFEM_updateJaco, CPFEM_dt)  ! parallel execution inside
           CPFEM_calc_done = .true.
         endif

!  translate from P and dP/dF to CS and dCS/dE
         Kirchhoff_bar = math_mul33x33(CPFEM_PK1_bar(:,:,CPFEM_in, cp_en),transpose(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en)))
         J_inverse  = 1.0_pReal/math_det3x3(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en))
         CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff_bar)
!
         H_bar = 0.0_pReal
         forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
           H_bar(i,j,k,l) = H_bar(i,j,k,l) + &
                          CPFEM_ffn1_bar(j,m,CPFEM_in,cp_en) * &
                          CPFEM_ffn1_bar(l,n,CPFEM_in,cp_en) * &
                          CPFEM_dPdF_bar(i,m,k,n,CPFEM_in,cp_en) - &
                          math_I3(j,l)*CPFEM_ffn1_bar(i,m,CPFEM_in,cp_en)*CPFEM_PK1_bar(k,m,CPFEM_in,cp_en) + &
                          0.5_pReal*(math_I3(i,k)*Kirchhoff_bar(j,l) + math_I3(j,l)*Kirchhoff_bar(i,k) + &
                                     math_I3(i,l)*Kirchhoff_bar(j,k) + math_I3(j,k)*Kirchhoff_bar(i,l))
         forall(i=1:3,j=1:3,k=1:3,l=1:3) &
            H_bar_sym(i,j,k,l)= 0.25_pReal*(H_bar(i,j,k,l) +H_bar(j,i,k,l) +H_bar(i,j,l,k) +H_bar(j,i,l,k))
         CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel3333to66(J_inverse*H_bar)
      endif
    case (3)    ! collect and return odd result
       CPFEM_Temperature(CPFEM_in,cp_en)  = Temperature
       CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)  = ffn
       CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) = ffn1
       CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_stress
       CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       CPFEM_calc_done = .false.

    case (4)    ! do nothing since we can recycle the former results (MARC specialty)
    case (5)    ! record consistent tangent at beginning of new increment (while recycling)
       CPFEM_jaco_knownGood = CPFEM_jaco_bar
    case (6)    ! restore consistent tangent after cutback
       CPFEM_jaco_bar = CPFEM_jaco_knownGood
 end select
!
! return the local stress and the jacobian from storage
 CPFEM_stress(1:CPFEM_ngens) = CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en)
 CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en)
!
 return
!
 END SUBROUTINE
!
!
!**********************************************************
!***      calculate the material point behaviour        ***
!**********************************************************
 SUBROUTINE CPFEM_MaterialPoint(&
     updateJaco,&     ! flag to initiate Jacobian updating
     CPFEM_dt)        ! Time increment (dt)
!
 use prec
 use debug
 use math, only: math_pDecomposition,math_RtoEuler,inDeg
 use IO,   only: IO_error
 use mesh, only: mesh_element, mesh_NcpElems, FE_Nips
 use material, only: homogenization_Ngrains,material_phase,material_volfrac
 use constitutive
 implicit none
!
 logical, intent(in) :: updateJaco
 real(pReal), intent(in) :: CPFEM_dt
 integer(pInt) g,i,e
 logical error
 real(pReal), dimension(3,3) :: U,R


!$OMP PARALLEL DO
 do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)       ! iterate over elements to be processed
   do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)     ! iterate over IPs of this element to be processed
     forall (g = 1:homogenization_Ngrains(mesh_element(3,e))) ! number of grains of this homogenization
          CPFEM_ffn(:,:,g,i,e)  = CPFEM_ffn_bar(:,:,i,e)      ! Taylor homogenization (why not using former ffn1??)
          CPFEM_ffn1(:,:,g,i,e) = CPFEM_ffn1_bar(:,:,i,e)     ! Taylor homogenization
     end forall
   enddo
 enddo
!$OMP END PARALLEL DO

 call SingleCrystallite(updateJaco,CPFEM_dt)

!******************************************************************************************************
! check convergence of homogenization if needed
!******************************************************************************************************

! calculate average quantities per ip and post results
!$OMP PARALLEL DO
 do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)    ! iterate over elements to be processed
   do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)  ! iterate over IPs of this element to be processed
     CPFEM_PK1_bar(:,:,i,e) = sum(CPFEM_PK1(:,:,:,i,e),3)/homogenization_Ngrains(mesh_element(3,e))
     if (updateJaco) &
       CPFEM_dPdF_bar(:,:,:,:,i,e) = &
         sum(CPFEM_dPdF(:,:,:,:,:,i,e),5)/homogenization_Ngrains(mesh_element(3,e))   ! add up crystallite stiffnesses (may have "holes" corresponding to former avg tangent)
     do g = 1,homogenization_Ngrains(mesh_element(3,e))
       call math_pDecomposition(CPFEM_Fe_new(:,:,g,i,e),U,R,error)      ! polar decomposition
       if (error) call IO_error(650,e,i,g)
       CPFEM_results(1,g,i,e) = material_phase(g,i,e)
       CPFEM_results(2,g,i,e) = material_volFrac(g,i,e)
       CPFEM_results(3:5,g,i,e) = math_RtoEuler(transpose(R))*inDeg     ! orientation
     enddo
   enddo
 enddo
!$OMP END PARALLEL DO

 return

 END SUBROUTINE


!********************************************************************
! Calculates the stress and jacobi (if wanted) for all or a single component
!********************************************************************
 subroutine SingleCrystallite(&
     updateJaco,& ! update of Jacobian required
     dt)          ! time increment

 use prec, only: pReal,pInt,pert_Fg,subStepMin, nCutback
 use debug
 use math
 use IO,   only: IO_error
 use mesh, only: mesh_element, FE_Nips
 use material, only: homogenization_Ngrains
 use constitutive

 implicit none

 character (len=128) msg
 logical updateJaco, allConverged
 real(preal) dt
 real(pReal), dimension(3,3) :: Fg_pert,Lp_pert, P_pert, Fp_pert, Fe_pert
 real(pReal), dimension(6) :: Tstar_v     
 real(pReal), dimension(constitutive_maxSizeState) :: state
 integer(pInt) g,i,e,k,l,iOuter,mySizeState

!$OMP PARALLEL DO
 do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)          ! iterate over elements to be processed
   do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)        ! iterate over IPs of this element to be processed
     forall (g = 1:homogenization_Ngrains(mesh_element(3,e)))    ! number of grains of this homogenization
       crystallite_converged(g,i,e) = .false.
       constitutive_state_new(g,i,e)%p = constitutive_state_old(g,i,e)%p
       CPFEM_Lp_new(:,:,g,i,e) = CPFEM_Lp_old(:,:,g,i,e)
     end forall
   end do
 end do
!$OMP END PARALLEL DO

 iOuter = 0_pInt
 allConverged = .false.
 
 do while (.not. allConverged)
   iOuter = iOuter + 1_pInt                                           ! count state integation loops
   if (iOuter > nOuter) call IO_error(600)                            ! too many loops required --> croak

!$OMP PARALLEL DO
   do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)             ! iterate over elements to be processed
     do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)           ! iterate over IPs of this element to be processed
       do g = 1,homogenization_Ngrains(mesh_element(3,e))             ! number of grains of this homogenization
         if (.not. crystallite_converged(g,i,e)) then
           call integrateStress(msg,CPFEM_Tstar_v(:,g,i,e),CPFEM_PK1(:,:,g,i,e), &
                                CPFEM_Fp_new(:,:,g,i,e),CPFEM_Fe_new(:,:,g,i,e),CPFEM_Lp_new(:,:,g,i,e), &
                                CPFEM_ffn1(:,:,g,i,e),dt,g,i,e)
           if (msg /= 'ok') call IO_error(610,e,i,g,msg)
         endif
       end do
     end do
   end do
!$OMP END PARALLEL DO

   allConverged = .true.                                              ! assume best case

!$OMP PARALLEL DO
   do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)             ! iterate over elements to be processed
     do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)           ! iterate over IPs of this element to be processed
       do g = 1,homogenization_Ngrains(mesh_element(3,e))             ! number of grains of this homogenization
         if (crystallite_converged(g,i,e)) cycle                      ! this one is already fine
         if (integrateState(CPFEM_Tstar_v(:,g,i,e),dt,g,i,e)) then    ! state integration now converged?
           crystallite_converged(g,i,e) = .true.
!$OMP CRITICAL (out)
           debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!$OMP END CRITICAL (out)
         else
           allConverged = .false.                                     ! this one requires additional round...
         endif
       end do
     end do
   end do
!$OMP END PARALLEL DO

 end do                                                               ! all crystallites converged
 
!$OMP PARALLEL DO
 do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)               ! iterate over elements to be processed
   do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)             ! iterate over IPs of this element to be processed
     forall (g = 1:homogenization_Ngrains(mesh_element(3,e))) &       ! number of grains of this homogenization
       CPFEM_results(CPFEM_Nresults+1:CPFEM_Nresults+constitutive_sizePostResults(g,i,e),g,i,e) = &
         constitutive_postResults(CPFEM_Tstar_v(:,g,i,e),CPFEM_Temperature(i,e),dt,g,i,e)
   end do
 end do
!$OMP END PARALLEL DO

 if(updateJaco) then                                                  ! Jacobian required

!$OMP CRITICAL (write2out)
   if (debugger) write (6,*) 'Jacobian calc'
!$OMP END CRITICAL (write2out)

!$OMP PARALLEL DO
   do e = CPFEM_execution_elem(1),CPFEM_execution_elem(2)             ! iterate over elements to be processed
     do i = CPFEM_execution_IP(1,e),CPFEM_execution_IP(2,e)           ! iterate over IPs of this element to be processed
       do g = 1,homogenization_Ngrains(mesh_element(3,e))             ! number of grains of this homogenization
         mySizeState = constitutive_sizeState(g,i,e)                  ! number of state variables for this grain
         state(1:mySizeState) = constitutive_state_new(g,i,e)%p       ! remember unperturbed, converged state
         do k = 1,3                                                   ! perturbation...
           do l = 1,3                                                 ! ...components
             Fg_pert = CPFEM_ffn1(:,:,g,i,e)                          ! initialize perturbed Fg
             Fg_pert(k,l) = Fg_pert(k,l) + pert_Fg                    ! perturb single component
             Lp_pert = CPFEM_Lp_new(:,:,g,i,e)                        ! initialize Lp
             Fp_pert = CPFEM_Fp_new(:,:,g,i,e)                        ! initialize Fp
             constitutive_state_new(g,i,e)%p = state(1:mySizeState)   ! initial guess from end of time step
             crystallite_converged(g,i,e) = .false.
             iOuter = 0_pInt
             do while(.not. crystallite_converged(g,i,e) .and. iOuter < nOuter)
               iOuter = iOuter + 1_pInt
               call integrateStress(msg,Tstar_v,P_pert,Fp_pert,Fe_pert,Lp_pert, Fg_pert,dt,g,i,e)
               if (msg /= 'ok') exit
               crystallite_converged(g,i,e) = integrateState(Tstar_v,dt,g,i,e)
             end do
             if (crystallite_converged(g,i,e)) &
               CPFEM_dPdF(:,:,k,l,g,i,e) = (P_pert-CPFEM_PK1(:,:,g,i,e))/pert_Fg ! constructing tangent dP_ij/dFg_kl only if valid forward difference
!$OMP CRITICAL (out)
             debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!$OMP END CRITICAL (out)
           end do
         end do
         constitutive_state_new(g,i,e)%p = state(1:mySizeState)      ! restore solution
       end do
     end do
   end do
!$OMP END PARALLEL DO
 endif

 return

 end subroutine


!********************************************************************
! Update the state for a single component
!********************************************************************
 function integrateState(&
   Tstar_v,&        ! stress
   dt,&             ! time increment
   g,&              ! grain number
   i,&              ! integration point number
   e&               ! element number
 )
 use prec, only: pReal,pInt,reltol_Outer
 use constitutive, only: constitutive_dotState,constitutive_sizeDotState,&
                         constitutive_state_old,constitutive_state_new
 
 logical integrateState

 integer(pInt) g,i,e,mySize
 real(pReal), dimension(6) :: Tstar_v
 real(pReal) dt
 real(pReal), dimension(constitutive_sizeDotState(g,i,e)) :: residuum

 mySize = constitutive_sizeDotState(g,i,e)
 residuum = constitutive_state_new(g,i,e)%p(1:mySize) - constitutive_state_old(g,i,e)%p(1:mySize) - &
            dt*constitutive_dotState(Tstar_v,CPFEM_Temperature(i,e),g,i,e)                        ! residuum from evolution of microstructure
 constitutive_state_new(g,i,e)%p(1:mySize) = constitutive_state_new(g,i,e)%p(1:mySize) - residuum ! update of microstructure
 integrateState = maxval(abs(residuum/constitutive_state_new(g,i,e)%p(1:mySize)),&
                             constitutive_state_new(g,i,e)%p(1:mySize) /= 0.0_pReal) < reltol_Outer
 return

 end function


!********************************************************************
! Calculates the stress for a single component
!********************************************************************
!***********************************************************************
!***     calculation of stress (P), stiffness (dPdF),                ***
!***     and announcement of any                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 subroutine integrateStress(&
     msg,&        ! return message
     Tstar_v,&    ! Stress vector
     P,&          ! first PK stress
     Fp_new,&     ! new plastic deformation gradient
     Fe_new,&     ! new "elastic" deformation gradient
     Lp,&         ! plastic velocity gradient
!
     Fg_new,&     ! new global deformation gradient
     dt,&         ! time increment
     g,&          ! grain number
     i,&          ! integration point number
     e)           ! element number

 use prec, only: pReal,pInt,pert_Fg,subStepMin, nCutback
 use debug
 use constitutive, only: constitutive_state_new
 use math
! use CPFEM
!
 implicit none
!
 character(len=*) msg
 logical error,success
 integer(pInt) e,i,g, nCutbacks, maxCutbacks
 real(pReal) Temperature
 real(pReal) dt,dt_aim,subFrac,subStep,det
 real(pReal), dimension(3,3)     :: Lp,Lp_interpolated,inv
 real(pReal), dimension(3,3)     :: Fg_current,Fg_new,Fg_aim,deltaFg
 real(pReal), dimension(3,3)     :: Fp_current,Fp_new
 real(pReal), dimension(3,3)     :: Fe_current,Fe_new
 real(pReal), dimension(3,3)     :: P
 real(pReal), dimension(6) :: Tstar_v

 deltaFg = Fg_new - CPFEM_ffn(:,:,g,i,e)
 subFrac = 0.0_pReal
 subStep = 1.0_pReal
 nCutbacks = 0_pInt
 maxCutbacks = 0_pInt
 Fg_current = CPFEM_ffn(:,:,g,i,e)                                 ! initialize to start of inc
 Fp_current = CPFEM_Fp_old(:,:,g,i,e)
 call math_invert3x3(Fp_current,inv,det,error)
 Fe_current = math_mul33x33(Fg_current,inv)

 success = .false.                                                 ! pretend cutback
 dt_aim = 0.0_pReal                                                ! prevent initial Lp interpolation
 Temperature = CPFEM_Temperature(i,e)

! begin the cutback loop
 do while (subStep > subStepMin)                                   ! continue until finished or too much cut backing
   if (success) then                                               ! wind forward
     Fg_current = Fg_aim
     Fe_current = Fe_new
     Fp_current = Fp_new
   elseif (dt_aim > 0.0_pReal) then
     call math_invert3x3(Fg_aim,inv,det,error)                     !  inv of Fg_aim
     Lp_interpolated = 0.5_pReal*Lp + &
                       0.5_pReal*(math_I3 -  math_mul33x33(Fp_current,&
                                  math_mul33x33(inv,Fe_current)))/dt_aim  ! interpolate Lp  and L
     if (debugger) then
!$OMP CRITICAL (write2out)
       write (6,*) 'Lp interpolation'
       write (6,'(a,/,3(3(f12.7,x)/))') 'from',Lp(1:3,:)
       write (6,'(a,/,3(3(f12.7,x)/))') 'to',Lp_interpolated(1:3,:)
!$OMP END CRITICAL (write2out)
     endif
     Lp = Lp_interpolated
   endif
!
   Fg_aim = Fg_current + subStep*deltaFg                           ! aim for Fg
   dt_aim = subStep*dt                                             ! aim for dt
   if (debugger) then
!$OMP CRITICAL (write2out)
     write (6,*) 'using these values'
     write (6,'(a,/,3(4(f9.3,x)/))') 'state new / MPa',constitutive_state_new(g,i,e)%p/1e6_pReal
     write (6,'(a,/,3(3(f12.7,x)/))') 'Fe current',Fe_current(1:3,:)
     write (6,'(a,/,3(3(f12.7,x)/))') 'Fp current',Fp_current(1:3,:)
     write (6,'(a,/,3(3(f12.7,x)/))') 'Lp (old=new guess)',Lp(1:3,:)
     write (6,'(a20,f,x,a2,x,f)') 'integrating from ',subFrac,'to',(subFrac+subStep)
!$OMP END CRITICAL (write2out)
   endif

   call TimeIntegration(msg,Lp,Fp_new,Fe_new,Tstar_v,P, Fg_aim,Fp_current,Temperature,dt_aim,g,i,e)

   if (msg == 'ok') then
     subFrac = subFrac + subStep
     subStep = min(1.0_pReal-subFrac, subStep*2.0_pReal) ! accelerate
     nCutbacks = 0_pInt                   ! reset cutback counter
     success = .true.                     ! keep current Lp
   else
     nCutbacks = nCutbacks + 1            ! record additional cutback
     maxCutbacks = max(nCutbacks,maxCutbacks) ! remember maximum number of cutbacks
     subStep = subStep / 2.0_pReal        ! cut time step in half
     success = .false.                    ! force Lp interpolation
   endif
 enddo  ! potential substepping
!
!$OMP CRITICAL (cutback)
 debug_cutbackDistribution(min(nCutback,maxCutbacks)+1) = debug_cutbackDistribution(min(nCutback,maxCutbacks)+1)+1
!$OMP END CRITICAL (cutback)

 return

 end subroutine

!
!***********************************************************************
!***     fully-implicit two-level time integration                   ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 SUBROUTINE TimeIntegration(&
     msg,&              ! return message
     Lpguess,&          ! guess of plastic velocity gradient
     Fp_new,&           ! new plastic deformation gradient
     Fe_new,&           ! new "elastic" deformation gradient
     Tstar_v,&          ! Stress vector
     P,&                ! 1st PK stress (taken as initial guess if /= 0)
     Fg_new,&           ! new total def gradient
     Fp_old,&           ! former plastic def gradient
     Temperature,&      ! temperature
     dt,&               ! time increment
     grain,&            ! grain number
     ip,&               ! integration point number
     cp_en &            ! element number
    )
    
 use prec
 use debug
 use mesh, only: mesh_element
 use constitutive, only: constitutive_microstructure,constitutive_homogenizedC,constitutive_LpAndItsTangent,&
                         constitutive_state_new
 use math
 use IO
 implicit none
!
 character(len=*) msg
 logical failed
 integer(pInt) cp_en, ip, grain
 integer(pInt) iInner,dummy, i,j,k,l,m,n
 real(pReal) dt, Temperature, det, p_hydro, leapfrog,maxleap
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(9,9) :: dLp,dTdLp,dRdLp,invdRdLp,eye2
 real(pReal), dimension(6,6) :: C_66
 real(pReal), dimension(3,3) :: Fg_new,Fp_new,invFp_new,Fp_old,invFp_old,Fe_new
 real(pReal), dimension(3,3) :: P
 real(pReal), dimension(3,3) :: Lp,Lpguess,Lpguess_old,Rinner,Rinner_old,A,B,BT,AB,BTA
 real(pReal), dimension(3,3,3,3) :: C

 msg = 'ok'  ! error-free so far
 eye2 = math_identity2nd(9)

 call math_invert3x3(Fp_old,invFp_old,det,failed) ! inversion of Fp_old
 if (failed) then
    msg = 'inversion Fp_old'
    return
 endif

 A = math_mul33x33(transpose(invFp_old), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_old)))

!$OMP CRITICAL (write2out)
 if (debugger) write (6,'(a,/,3(3(f12.7,x)/))') 'Fg to be calculated',Fg_new
!$OMP END CRITICAL (write2out)

 call constitutive_microstructure(Temperature,grain,ip,cp_en)
 C_66 = constitutive_homogenizedC(grain,ip,cp_en)
 C = math_Mandel66to3333(C_66)       ! 4th rank elasticity tensor

 iInner = 0_pInt
 leapfrog = 1.0_pReal                ! correction as suggested by invdRdLp-step
 maxleap = 1024.0_pReal              ! preassign maximum acceleration level

 Lpguess_old = Lpguess               ! consider present Lpguess good

Inner: do              ! inner iteration: Lp
   iInner = iInner+1
   if (iInner > nInner) then           ! too many loops required
     Lpguess = Lpguess_old             ! do not trust the last update but resort to former one
     msg = 'limit Inner iteration'
     return
   endif

   B = math_i3 - dt*Lpguess
   BT = transpose(B)
   AB = math_mul33x33(A,B)
   BTA = math_mul33x33(BT,A)
   Tstar_v = 0.5_pReal*math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB)-math_I3))
   p_hydro=(Tstar_v(1)+Tstar_v(2)+Tstar_v(3))/3.0_pReal
   forall(i=1:3) Tstar_v(i) = Tstar_v(i)-p_hydro                ! subtract hydrostatic pressure
   call constitutive_LpAndItsTangent(Lp,dLp, Tstar_v,Temperature,grain,ip,cp_en)
   Rinner = Lpguess - Lp                                        ! update current residuum

   if (.not.(any(Rinner/=Rinner)) .and. &                       ! exclude any NaN in residuum
       ( ( maxval(abs(Rinner)) < abstol_Inner) .or. &            ! below abs tol .or.
         ( any(abs(dt*Lpguess) > relevantStrain) .and. &        ! worth checking? .and.
           maxval(abs(Rinner/Lpguess),abs(dt*Lpguess) > relevantStrain) < reltol_Inner &  ! below rel tol
         ) &
       )   &
      )    &
     exit Inner                                                ! convergence
!
!          check for acceleration/deceleration in Newton--Raphson correction
!
   if (any(Rinner/=Rinner) .and. &                              ! NaN occured at regular speed
       leapfrog == 1.0) then
     Lpguess = Lpguess_old                                     ! restore known good guess
     msg = 'NaN present'                                       ! croak for cutback
     return

   elseif (leapfrog > 1.0_pReal .and. &                         ! at fast pace ?
           (sum(Rinner*Rinner) > sum(Rinner_old*Rinner_old) .or. &  ! worse residuum
            sum(Rinner*Rinner_old) < 0.0_pReal) .or. &              ! residuum changed sign (overshoot)
           any(Rinner/=Rinner) ) then                              ! NaN
     maxleap = 0.5_pReal * leapfrog                             ! limit next acceleration
     leapfrog = 1.0_pReal                                       ! grinding halt

   else                                                         ! better residuum
     dTdLp = 0.0_pReal                                          ! calc dT/dLp
     forall (i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
       dTdLp(3*(i-1)+j,3*(k-1)+l) = dTdLp(3*(i-1)+j,3*(k-1)+l) + &
                                    C(i,j,l,n)*AB(k,n)+C(i,j,m,l)*BTA(m,k)
     dTdLp = -0.5_pReal*dt*dTdLp
     dRdLp = eye2 - math_mul99x99(dLp,dTdLp)                    ! calc dR/dLp
     invdRdLp = 0.0_pReal
     call math_invert(9,dRdLp,invdRdLp,dummy,failed)            ! invert dR/dLp --> dLp/dR
     if (failed) then
       msg = 'inversion dR/dLp'
       if (debugger) then
!$OMP CRITICAL (write2out)
                 write (6,*) msg
                 write (6,'(a,/,9(9(e9.3,x)/))') 'dRdLp', dRdLp(1:9,:)
                 write (6,'(a,/,3(4(f9.3,x)/))') 'state_new / MPa',constitutive_state_new(grain,ip,cp_en)%p/1e6_pReal
                 write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
                 write (6,'(a,/,3(3(e12.7,x)/))') 'Lp',Lp(1:3,:)
                 write (6,'(a,/,6(f9.3,x))') 'Tstar / MPa',Tstar_v/1e6_pReal
!$OMP END CRITICAL (write2out)
       endif 
       return
     endif
!
     Rinner_old = Rinner                                        ! remember current residuum
     Lpguess_old = Lpguess                                      ! remember current Lp guess 
     if (iInner > 1 .and. leapfrog < maxleap) leapfrog = 2.0_pReal * leapfrog   ! accelerate if ok
   endif
!
   Lpguess = Lpguess_old                                        ! start from current guess                                       
   Rinner  = Rinner_old                                         ! use current residuum
   forall (i=1:3,j=1:3,k=1:3,l=1:3) &                           ! leapfrog to updated Lpguess 
     Lpguess(i,j) = Lpguess(i,j) - leapfrog*invdRdLp(3*(i-1)+j,3*(k-1)+l)*Rinner(k,l)
 enddo Inner
!
!$OMP CRITICAL (in)
 debug_InnerLoopDistribution(iInner) = debug_InnerLoopDistribution(iInner)+1
!$OMP END CRITICAL (in)
 invFp_new = math_mul33x33(invFp_old,B)
 call math_invert3x3(invFp_new,Fp_new,det,failed)
 if (failed) then
   msg = 'inversion Fp_new^-1'
   return
 endif

 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal)       ! regularize Fp by det = det(InvFp_new) !!
 forall (i=1:3) Tstar_v(i) = Tstar_v(i) + p_hydro ! add hydrostatic component back
 Fe_new = math_mul33x33(Fg_new,invFp_new)              ! calc resulting Fe
 P = math_mul33x33(Fe_new,math_mul33x33(math_Mandel6to33(Tstar_v),transpose(invFp_new)))    ! first PK stress

 return
!
 END SUBROUTINE
!
 END MODULE
!##############################################################