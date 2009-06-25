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
 real(pReal), dimension (:,:),        allocatable :: CPFEM_Temperature
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_ffn_bar
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_ffn1_bar
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_PK1_bar
 real(pReal), dimension (:,:,:,:,:,:),allocatable :: CPFEM_dPdF_bar
 real(pReal), dimension (:,:,:),      allocatable :: CPFEM_stress_bar
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_jaco_bar
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_jaco_knownGood
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_results
 real(pReal), dimension (:,:,:,:,:),  allocatable :: CPFEM_Fp_old
 real(pReal), dimension (:,:,:,:,:),  allocatable :: CPFEM_Fp_new
 real(pReal), parameter :: CPFEM_odd_stress = 1e15_pReal, CPFEM_odd_jacobian = 1e50_pReal
 integer(pInt) :: CPFEM_Nresults   = 4_pInt    ! three Euler angles plus volume fraction
 logical :: CPFEM_init_done = .false.          ! remember if init has been done already
 logical :: CPFEM_calc_done = .false.          ! remember if first IP has already calced the results
!
 real(pReal), dimension (:,:,:,:),    allocatable :: GIA_rVect_new    ! boundary relaxation vectors
 real(pReal), dimension (:,:,:,:),    allocatable :: GIA_rVect_old    ! boundary relaxation vectors
 real(pReal), dimension (:,:),        allocatable :: GIA_bNorm    ! grain boundary normals
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
 use mesh
 use constitutive
!
 implicit none
!
 real(pReal) Temperature
 integer(pInt) e,i,g,b
!
!    *** mpie.marc parameters ***
 allocate(CPFEM_Temperature   (mesh_maxNips,mesh_NcpElems)) ; CPFEM_Temperature = Temperature
 allocate(CPFEM_ffn_bar   (3,3,mesh_maxNips,mesh_NcpElems))
 forall(e=1:mesh_NcpElems,i=1:mesh_maxNips) CPFEM_ffn_bar(:,:,i,e)             = math_I3
 allocate(CPFEM_ffn1_bar  (3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1_bar   = CPFEM_ffn_bar
 allocate(CPFEM_PK1_bar   (3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_PK1_bar    = 0.0_pReal
 allocate(CPFEM_dPdF_bar(3,3,3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF_bar = 0.0_pReal
 allocate(CPFEM_stress_bar(6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_stress_bar = 0.0_pReal
 allocate(CPFEM_jaco_bar(6,6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_jaco_bar   = 0.0_pReal
 allocate(CPFEM_jaco_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_jaco_knownGood = 0.0_pReal
!
!    *** User defined results !!! MISSING incorporate consti_Nresults ***
 allocate(CPFEM_results(CPFEM_Nresults+constitutive_maxNresults,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_results = 0.0_pReal
!
!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fp_new = 0.0_pReal
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall (e=1:mesh_NcpElems,i=1:mesh_maxNips,g=1:constitutive_maxNgrains) &
   CPFEM_Fp_old(:,:,g,i,e) = math_EulerToR(constitutive_EulerAngles(:,g,i,e))  ! plastic def gradient reflects init orientation
!
 allocate(GIA_rVect_new(3,12,mesh_maxNips,mesh_NcpElems)) ; GIA_rVect_new = 0.0_pReal
 allocate(GIA_rVect_old(3,12,mesh_maxNips,mesh_NcpElems)) ; GIA_rVect_old = 0.0_pReal
 allocate(GIA_bNorm(3,12))                                ; GIA_bNorm = 0.0_pReal
 do b = 1,4 
   GIA_bNorm(1,b)   = 1.0_pReal
   GIA_bNorm(2,b+4) = 1.0_pReal
   GIA_bNorm(3,b+8) = 1.0_pReal
 enddo
!
!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)

 write(6,*)
 write(6,*) 'CPFEM Initialization'
 write(6,*)
 write(6,*) 'CPFEM_Temperature:   ', shape(CPFEM_Temperature)
 write(6,*) 'CPFEM_ffn_bar:       ', shape(CPFEM_ffn_bar)
 write(6,*) 'CPFEM_ffn1_bar:      ', shape(CPFEM_ffn1_bar)
 write(6,*) 'CPFEM_PK1_bar:       ', shape(CPFEM_PK1_bar)
 write(6,*) 'CPFEM_dPdF_bar:      ', shape(CPFEM_dPdF_bar)
 write(6,*) 'CPFEM_stress_bar:    ', shape(CPFEM_stress_bar)
 write(6,*) 'CPFEM_jaco_bar:      ', shape(CPFEM_jaco_bar)
 write(6,*) 'CPFEM_jaco_knownGood: ', shape(CPFEM_jaco_knownGood)
 write(6,*) 'CPFEM_results:       ', shape(CPFEM_results)
 write(6,*) 'CPFEM_Fp_old:        ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:        ', shape(CPFEM_Fp_new)
!
 write(6,*) 'GIA_rVect_new:       ', shape(GIA_rVect_new)
 write(6,*) 'GIA_rVect_old:       ', shape(GIA_rVect_old)
 write(6,*) 'GIA_bNorm:           ', shape(GIA_bNorm)
 write(6,*)
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
 use mesh, only: mesh_init,mesh_FEasCP, mesh_NcpElems, FE_Nips, FE_mapElemtype, mesh_element
 use lattice, only: lattice_init
 use constitutive, only: constitutive_init,constitutive_state_old,constitutive_state_new,material_Cslip_66
 implicit none
!
 integer(pInt) CPFEM_en, CPFEM_in, cp_en, CPFEM_ngens, i,j,k,l,m,n, e
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
   call mesh_init()
   call lattice_init()
   call constitutive_init()
   call CPFEM_init(Temperature)			
   CPFEM_init_done = .true.
 endif
!
 cp_en = mesh_FEasCP('elem',CPFEM_en)
  if (cp_en == 1 .and. CPFEM_in == 1) then
!$OMP CRITICAL (write2out)
    write(6,'(a6,x,i4,x,a4,x,i4,x,a10,x,f8.4,x,a10,x,i2,x,a10,x,i2,x,a10,x,i2,x,a10,x,i2)') &
    'elem',cp_en,'IP',CPFEM_in,&
    'theTime',theTime,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,&
    'mode',CPFEM_mode
!$OMP END CRITICAL (write2out)
  endif
!
 select case (CPFEM_mode)
    case (2,1)     ! regular computation (with aging of results)
       if (.not. CPFEM_calc_done) then                ! puuh, me needs doing all the work...
!$OMP CRITICAL (write2out)
           write (6,*) 'puuh me needs doing all the work', cp_en
!$OMP END CRITICAL (write2out)
           if (CPFEM_mode == 1) then                  ! age results at start of new increment
             CPFEM_Fp_old            = CPFEM_Fp_new
             constitutive_state_old  = constitutive_state_new
             GIA_rVect_old           = GIA_rVect_new
!$OMP CRITICAL (write2out)
             write (6,*) '#### aged results'
!$OMP END CRITICAL (write2out)
           endif
           debug_cutbackDistribution = 0_pInt         ! initialize debugging data
           debug_InnerLoopDistribution = 0_pInt
           debug_OuterLoopDistribution = 0_pInt
!
           do e=1,mesh_NcpElems                       ! ## this shall be done in a parallel loop in the future ##
               do i=1,FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
                  debugger = (e==1 .and. i==1)        ! switch on debugging for first IP in first element
                  call CPFEM_MaterialPoint(CPFEM_updateJaco, CPFEM_dt, i, e)
               enddo
           enddo
           call debug_info()                          ! output of debugging/performance statistics
           CPFEM_calc_done = .true.                   ! now calc is done
         endif    
!       translate from P and dP/dF to CS and dCS/dE
!!$OMP CRITICAL (evilmatmul)
       Kirchhoff_bar = math_mul33x33(CPFEM_PK1_bar(:,:,CPFEM_in, cp_en),transpose(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en)))
!!$OMP END CRITICAL (evilmatmul)
       J_inverse  = 1.0_pReal/math_det3x3(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en))
       CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff_bar)
!
       H_bar = 0.0_pReal
       forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
         H_bar(i,j,k,l) = H_bar(i,j,k,l) + &
                          (CPFEM_ffn1_bar(j,m,CPFEM_in,cp_en)*CPFEM_ffn1_bar(l,n,CPFEM_in,cp_en)*CPFEM_dPdF_bar(i,m,k,n,CPFEM_in,cp_en) - &
                           math_I3(j,l)*CPFEM_ffn1_bar(i,m,CPFEM_in,cp_en)*CPFEM_PK1_bar(k,m,CPFEM_in,cp_en)) + &
                          0.5_pReal*(math_I3(i,k)*Kirchhoff_bar(j,l) + math_I3(j,l)*Kirchhoff_bar(i,k) + &
                                     math_I3(i,l)*Kirchhoff_bar(j,k) + math_I3(j,k)*Kirchhoff_bar(i,l))
       forall(i=1:3,j=1:3,k=1:3,l=1:3) &
          H_bar_sym(i,j,k,l)= 0.25_pReal*(H_bar(i,j,k,l) +H_bar(j,i,k,l) +H_bar(i,j,l,k) +H_bar(j,i,l,k))
       CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel3333to66(J_inverse*H_bar)
!
    case (3)    ! collect and return odd result
       CPFEM_Temperature(CPFEM_in,cp_en)  = Temperature
       CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)  = ffn
       CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) = ffn1
       CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_stress
       CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       CPFEM_calc_done = .false.
    case (4)    ! do nothing since we can recycle the former results (MARC specialty)
    case (5)    ! record consistent tangent at beginning of new increment
       CPFEM_jaco_knownGood = CPFEM_jaco_bar
    case (6)    ! restore consistent tangent after cutback
       CPFEM_jaco_bar = CPFEM_jaco_knownGood
 end select
!
! return the local stress and the jacobian from storage
 CPFEM_stress(1:CPFEM_ngens) = CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en)
 CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en)
! if (cp_en == 1 .and. CPFEM_in == 1) write (6,*) 'stress',CPFEM_stress
! if (cp_en == 1 .and. CPFEM_in == 1 .and. CPFEM_updateJaco) write (6,*) 'stiffness',CPFEM_jaco
! if (cp_en == 1 .and. CPFEM_in == 1) write (6,*) 'vector',GIA_rVect_new(:,:,1,1)
! 
 return
!
 END SUBROUTINE
!
!**********************************************************
!***      calculate the material point behaviour        ***
!**********************************************************
 SUBROUTINE CPFEM_MaterialPoint(&
     updateJaco,&     ! flag to initiate Jacobian updating
     CPFEM_dt,&       ! Time increment (dt)
     CPFEM_in,&       ! Integration point number
     cp_en)           ! Element number
!
 use prec
 use FEsolving, only: theCycle
 use debug
 use math, only: math_pDecomposition,math_RtoEuler,inDeg,math_I3,math_invert3x3,math_permut,math_invert,math_delta
 use IO,   only: IO_error
 use mesh, only: mesh_element
 use crystallite
 use constitutive
 implicit none
!
 character(len=128) msg
 integer(pInt) cp_en,CPFEM_in,grain,max_cutbacks,i,j,k,l,m,n,iBoun,NRiter,dummy,ii,jj,kk,ll,ip,jp
 logical updateJaco,error,NRconvergent,failed
 real(pReal) CPFEM_dt,volfrac,dTime,shMod,C_kb,resNorm,resMax,subStep,subFrac,temp1,temp2
 real(pReal), dimension(3,3)     :: F0_bar,F1_bar,dF_bar,PK1_per,F1_per
 real(pReal), dimension(3,3)     :: U,R
 real(pReal), dimension(3,3,8)      :: PK1,Fp0,Fp1,Fe1,F1,F0
 real(pReal), dimension(3,3,12)     :: GPK1,GF1,Nye,GRB1
 real(pReal), dimension(3,3,3,3,8)  :: dPdF
 real(pReal), dimension(3,3,3,3,12) :: dRdX1
 real(pReal), dimension(36)         :: var,res
 real(pReal), dimension(36,36)      :: dresdvar,dvardres
 real(pReal), dimension(3,12)       :: rx,rVect
 real(pReal), dimension(12)         :: NyeNorm
 real(pReal), dimension(constitutive_maxNstatevars,8) :: state0,state1
!
 if (texture_Ngrains(mesh_element(4,cp_en)) /= 8_pInt) then
   call IO_error(800)
   return
 endif
!
 CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) = 0.0_pReal                       ! zero out average first PK stress
 if (updateJaco) CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) = 0.0_pReal  ! zero out average consistent tangent
!
! ------------- GIA loop --------------------
!
! collect information
 shMod  = 0.2_pReal*(material_C11(1) - material_C12(1)) + 0.3_pReal*material_C44(1) ! equivalent shear modulus
 C_kb   = material_bg(1)*shMod/material_GrainSize(1)          ! equivalent boundary stiffness
!
 F0_bar = CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)                   ! effective deformation gradient at t_n
 state0 = constitutive_state_old(:,:,CPFEM_in,cp_en)          ! state variables at t_n
 Fp0    = CPFEM_Fp_old(:,:,:,CPFEM_in,cp_en)                  ! grain plastic def. gradient at t_n
 rVect  = GIA_rVect_old(:,:,CPFEM_in,cp_en)                   ! relaxation vectors from previous convergent step
!
 dF_bar = CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) - CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)   ! deformation gradient increment
 subFrac = 0.0_pReal
 subStep = 1.0_pReal
!
! Substepping procedure to improve N-R iteration
 SubStepping: do
 dTime  = subStep*CPFEM_dt
 call GIA_RelaxedDeformation(F0,F0_bar,rVect)                 ! def. gradient of indiv. grains at t_n
 F1_bar = F0_bar + subStep*dF_bar                             ! effective def. gradient at t_n+1
 forall (iBoun=1:12,i=1:3) var(3_pInt*(iBoun-1_pInt)+i) = rVect(i,iBoun)  ! primary variable: relaxation vector
!
! Newton-Raphson iteration block
 NRiter = 1_pInt
 NRIteration: do
   forall (iBoun=1:12,i=1:3) rx(i,iBoun) = var(3_pInt*(iBoun-1_pInt)+i)   ! relaxation vectors (guess)
!
!   deformation gradients of grains at t_n+1 (guess)
   call GIA_RelaxedDeformation(F1,F1_bar,rx)
!
! -------------- grain loop -----------------
   do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
     call SingleCrystallite(msg,PK1(:,:,grain),dPdF(:,:,:,:,grain),&
                      CPFEM_results(CPFEM_Nresults+1:CPFEM_Nresults+constitutive_Nresults(grain,CPFEM_in,cp_en),&
                                    grain,CPFEM_in,cp_en),&
                      Fp1(:,:,grain),Fe1(:,:,grain),state1(:,grain),&   ! output up to here
                      dTime,cp_en,CPFEM_in,grain,.true.,&
                      CPFEM_Temperature(CPFEM_in,cp_en),F1(:,:,grain),F0(:,:,grain),Fp0(:,:,grain),state0(:,grain))
     if (msg /= 'ok') then                    ! solution not reached --> exit NRIteration
!$OMP CRITICAL (write2out)
       write(6,*) 'GIA: grain loop failed to converge @ EL:',cp_en,' IP:',CPFEM_in
!$OMP END CRITICAL (write2out)
       NRconvergent = .false.
       exit NRiteration
     endif
   enddo    ! grain loop
!
!   calculate the deformation jump and stress jump across the boundaries
   call GIA_BoundaryJump(GF1,F1)
   call GIA_BoundaryJump(GPK1,PK1)
!
!   compute the Nye tensor at the boundary
   Nye     = 0.0_pReal
   NyeNorm = 0.0_pReal
   do iBoun = 1,12
     do i = 1,3
     do j = 1,3
       do k = 1,3
       do l = 1,3
         Nye(i,j,iBoun) = Nye(i,j,iBoun) - 0.5_pReal*math_permut(j,k,l)*GIA_bNorm(k,iBoun)*GF1(i,l,iBoun)
       enddo
       enddo
       NyeNorm(iBoun) = NyeNorm(iBoun) + Nye(i,j,iBoun)*Nye(i,j,iBoun)
     enddo
     enddo
     NyeNorm(iBoun) = sqrt(NyeNorm(iBoun))
     if (NyeNorm(iBoun) > 1.0e-8_pReal) Nye(:,:,iBoun) = Nye(:,:,iBoun)/NyeNorm(iBoun)
   enddo
!
!   compute the stress-like penalty at the boundary
   GRB1 = 0.0_pReal
   do iBoun = 1,12
     do i = 1,3
     do j = 1,3
       do k = 1,3
       do l = 1,3
         GRB1(i,j,iBoun) = GRB1(i,j,iBoun) + Nye(i,k,iBoun)*GIA_bNorm(l,iBoun)*math_permut(k,l,j)
       enddo
       enddo
     enddo
     enddo
     GRB1(:,:,iBoun) = 0.5_pReal*(C_kb + C_kb)*GRB1(:,:,iBoun)
   enddo
!
!   compute the resiudal of stress at the boundary
   res     = 0.0_pReal
   resNorm = 0.0_pReal
   do iBoun = 1,12
     do j = 1,3
       do i = 1,3
         res(3_pInt*(iBoun-1_pInt)+j) = res(3_pInt*(iBoun-1_pInt)+j) - &
                                        GIA_bNorm(i,iBoun)*(GPK1(i,j,iBoun) - GRB1(i,j,iBoun))
       enddo
       resNorm = resNorm + res(3_pInt*(iBoun-1_pInt)+j)*res(3_pInt*(iBoun-1_pInt)+j)
     enddo
   enddo
   resNorm = sqrt(resNorm)
!
   if (debugger) then
!$OMP CRITICAL (write2out)
    write(6,'(x,a,i3,a,i3,a,i3,a,e10.4)')'EL:',cp_en,' IP:',CPFEM_in,' Iter:',NRiter,' RNorm:',resNorm
!$OMP END CRITICAL (write2out)
   if (NRiter == 1_pInt) resMax = resNorm
   if ((resNorm < resToler*resMax) .or. (resNorm < resAbsol)) then     ! resNorm < tolerance ===> convergent
     NRconvergent = .true.
     exit NRiteration
   elseif ((NRiter > NRiterMax) .or. (resNorm > resBound*resMax)) then ! resNorm > up. bound ===> substepping
     NRconvergent = .false.
     exit NRiteration
   else                                                                ! update the residual
     dRdX1 = 0.0_pReal
     do iBoun = 1,12
       if (NyeNorm(iBoun) < 1.0e-8_pReal) NyeNorm(iBoun) = 1.0e-8_pReal
       do i = 1,3
       do j = 1,3
       do k = 1,3
       do l = 1,3
         temp1 = 0.0_pReal
         temp2 = 0.0_pReal
         do ii = 1,3
         do jj = 1,3
         do kk = 1,3
           temp1 = temp1 + GIA_bNorm(jj,iBoun)*math_permut(ii,jj,j)*math_delta(i,k)* &
                           GIA_bNorm(kk,iBoun)*math_permut(ii,kk,l)
           do ll = 1,3
             temp2 = temp2 + Nye(i,ii,iBoun)*GIA_bNorm(jj,iBoun)*math_permut(ii,jj,j)* &
                             Nye(k,kk,iBoun)*GIA_bNorm(ll,iBoun)*math_permut(kk,ll,l)
           enddo
         enddo
         enddo
         enddo
         dRdX1(i,j,k,l,iBoun) = 0.25_pReal*(C_kb + C_kb)*(temp1 - temp2)/NyeNorm(iBoun)
       enddo
       enddo
       enddo
       enddo
     enddo
     call GIA_JacobianMatrix(dresdvar,dPdF,dRdX1)
     dvardres = 0.0_pReal
     call math_invert(36,dresdvar,dvardres,dummy,failed)
     if (failed) then
!$OMP CRITICAL (write2out)
       write(6,*) 'GIA: failed to invert the Jacobian @ EL:',cp_en,' IP:',CPFEM_in
!$OMP END CRITICAL (write2out)
       NRconvergent = .false.
       exit NRiteration
     endif
     forall (i=1:36,j=1:36) var(i) = var(i) - dvardres(i,j)*res(j)
   endif
!
   NRiter = NRiter + 1_pInt
 enddo NRIteration ! End of N-R iteration blok
!
 if (.not. NRconvergent) then
   subStep = 0.5_pReal*subStep
 else
   subFrac = subFrac + subStep
   subStep = 1.0_pReal - subFrac
   Fp0    = Fp1
   F0_bar = F1_bar
   state0 = state1
   rVect  = rx
 endif
!
 if (subStep < subStepMin) exit SubStepping
 enddo SubStepping ! End of substepping blok
!
! ------------- GIA loop (end) --------------
!
! return to the general subroutine when convergence is not reached
 if (.not. NRconvergent) then
!$OMP CRITICAL (write2out)
   write(6,'(x,a)') 'GIA: convergence is not reached @ EL:',cp_en,' IP:',CPFEM_in
!$OMP END CRITICAL (write2out)
   call IO_error(600)
   return
 endif
!
! updates all variables, deformation gradients, and vectors
 GIA_rVect_new(:,:,CPFEM_in,cp_en)          = rVect
 CPFEM_Fp_new(:,:,:,CPFEM_in,cp_en)         = Fp1
 constitutive_state_new(:,:,CPFEM_in,cp_en) = state1
!
!   compute the effective stress and consistent tangent
 do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
   volfrac = constitutive_matVolFrac(grain,CPFEM_in,cp_en)*constitutive_texVolFrac(grain,CPFEM_in,cp_en)
   CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) = CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) + &
                                       volfrac*PK1(:,:,grain)                  ! average Cauchy stress
!
!   update results plotted in MENTAT
   call math_pDecomposition(Fe1(:,:,grain),U,R,error) ! polar decomposition
   if (error) then
!$OMP CRITICAL (write2out)
     write(6,*) Fe1(:,:,grain)
     write(6,*) 'polar decomposition'
     write(6,*) 'Grain:             ',grain
     write(6,*) 'Integration point: ',CPFEM_in
     write(6,*) 'Element:           ',mesh_element(1,cp_en)
!$OMP END CRITICAL (write2out)
     call IO_error(650)
     return
   endif
   CPFEM_results(1:3,grain,CPFEM_in,cp_en) = math_RtoEuler(transpose(R))*inDeg ! orientation
   CPFEM_results(4  ,grain,CPFEM_in,cp_en) = volfrac                           ! volume fraction of orientation
 enddo
!
 if (theCycle >= 0_pInt) then
   forall (grain=1:texture_Ngrains(mesh_element(4,cp_en))) &
      CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) = CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) + volfrac*dPdF(:,:,:,:,grain)
 else
   do ip = 1,3
   do jp = 1,3
     F1_per = F1_bar
     F1_per(ip,jp) = F1_per(ip,jp) + 1.0e-5_pReal
     forall (iBoun=1:12,i=1:3) var(3_pInt*(iBoun-1_pInt)+i) = rVect(i,iBoun)
     NRiter = 1_pInt
!
     NRPerturbation: do
       forall (iBoun=1:12,i=1:3) rx(i,iBoun) = var(3_pInt*(iBoun-1_pInt)+i)   ! relaxation vectors (guess)
       call GIA_RelaxedDeformation(F1,F1_bar,rx)
       do grain = 1,8
         call SingleCrystallite(msg,PK1(:,:,grain),dPdF(:,:,:,:,grain),&
                      CPFEM_results(CPFEM_Nresults+1:CPFEM_Nresults+constitutive_Nresults(grain,CPFEM_in,cp_en),&
                                    grain,CPFEM_in,cp_en),&
                      Fp1(:,:,grain),Fe1(:,:,grain),state1(:,grain),&   ! output up to here
                      dTime,cp_en,CPFEM_in,grain,.true.,&
                      CPFEM_Temperature(CPFEM_in,cp_en),F1(:,:,grain),F0(:,:,grain),Fp0(:,:,grain),state0(:,grain))
         if (msg /= 'ok') then                    ! solution not reached --> exit NRIteration
!$OMP CRITICAL (write2out)
           write(6,*) 'GIA: perturbation grain loop failed to converge within allowable step-size'
!$OMP END CRITICAL (write2out)
           NRconvergent = .false.
           exit NRPerturbation
         endif
       enddo
       call GIA_BoundaryJump(GF1,F1)
       call GIA_BoundaryJump(GPK1,PK1)
!
       Nye     = 0.0_pReal
       NyeNorm = 0.0_pReal
       do iBoun = 1,12
         do i = 1,3
         do j = 1,3
           do k = 1,3
           do l = 1,3
             Nye(i,j,iBoun) = Nye(i,j,iBoun) - 0.5_pReal*math_permut(j,k,l)*GIA_bNorm(k,iBoun)*GF1(i,l,iBoun)
           enddo
           enddo
           NyeNorm(iBoun) = NyeNorm(iBoun) + Nye(i,j,iBoun)*Nye(i,j,iBoun)
         enddo
         enddo
         NyeNorm(iBoun) = sqrt(NyeNorm(iBoun))
         if (NyeNorm(iBoun) > 1.0e-8_pReal) Nye(:,:,iBoun) = Nye(:,:,iBoun)/NyeNorm(iBoun)
       enddo
!
       GRB1 = 0.0_pReal
       do iBoun = 1,12
         do i = 1,3
         do j = 1,3
           do k = 1,3
           do l = 1,3
             GRB1(i,j,iBoun) = GRB1(i,j,iBoun) + Nye(i,k,iBoun)*GIA_bNorm(l,iBoun)*math_permut(k,l,j)
           enddo
           enddo
         enddo
         enddo
         GRB1(:,:,iBoun) = 0.5_pReal*(C_kb + C_kb)*GRB1(:,:,iBoun)
       enddo
!
       res     = 0.0_pReal
       resNorm = 0.0_pReal
       do iBoun = 1,12
         do j = 1,3
           do i = 1,3
             res(3_pInt*(iBoun-1_pInt)+j) = res(3_pInt*(iBoun-1_pInt)+j) - &
                                            GIA_bNorm(i,iBoun)*(GPK1(i,j,iBoun) - GRB1(i,j,iBoun))
           enddo
           resNorm = resNorm + res(3_pInt*(iBoun-1_pInt)+j)*res(3_pInt*(iBoun-1_pInt)+j)
         enddo
       enddo
       resNorm = sqrt(resNorm)
!
!       if (debugger) then
!!$OMP CRITICAL (write2out)
!            write(6,'(x,a,i3,a,i3,a,i3,a,i3,a,e10.4)')'EL = ',cp_en,':IP = ',CPFEM_in,':pert = ',3*(ip-1)+jp,':Iter = ',NRiter,':RNorm = ',resNorm
!!$OMP END CRITICAL (write2out)
!       endif
       if (NRiter == 1_pInt) resMax = resNorm
       if ((resNorm < resToler*resMax) .or. (resNorm < resAbsol)) then     ! resNorm < tolerance ===> convergent
         NRconvergent = .true.
         exit NRPerturbation
       elseif ((NRiter > NRiterMax) .or. (resNorm > resBound*resMax)) then ! resNorm > up. bound ===> substepping
         NRconvergent = .false.
         exit NRPerturbation
       else                                                                ! update the residual
         dRdX1 = 0.0_pReal
         do iBoun = 1,12
           if (NyeNorm(iBoun) < 1.0e-8_pReal) NyeNorm(iBoun) = 1.0e-8_pReal
           do i = 1,3
           do j = 1,3
           do k = 1,3
           do l = 1,3
             temp1 = 0.0_pReal
             temp2 = 0.0_pReal
             do ii = 1,3
             do jj = 1,3
             do kk = 1,3
               temp1 = temp1 + GIA_bNorm(jj,iBoun)*math_permut(ii,jj,j)*math_delta(i,k)* &
                               GIA_bNorm(kk,iBoun)*math_permut(ii,kk,l)
               do ll = 1,3
                 temp2 = temp2 + Nye(i,ii,iBoun)*GIA_bNorm(jj,iBoun)*math_permut(ii,jj,j)* &
                                 Nye(k,kk,iBoun)*GIA_bNorm(ll,iBoun)*math_permut(kk,ll,l)
               enddo
             enddo
             enddo
             enddo
             dRdX1(i,j,k,l,iBoun) = 0.25_pReal*(C_kb + C_kb)*(temp1 - temp2)/NyeNorm(iBoun)
           enddo
           enddo
           enddo
           enddo
         enddo
         call GIA_JacobianMatrix(dresdvar,dPdF,dRdX1)
         dvardres = 0.0_pReal
         call math_invert(36,dresdvar,dvardres,dummy,failed)
         if (failed) then
!$OMP CRITICAL (write2out)
           write(6,*) 'GIA: perturbation failed to invert the Jacobian'
!$OMP END CRITICAL (write2out)
           NRconvergent = .false.
           exit NRPerturbation
         endif
         forall (i=1:36,j=1:36) var(i) = var(i) - dvardres(i,j)*res(j)
       endif
       NRiter = NRiter + 1_pInt
     enddo NRPerturbation ! End of N-R iteration blok
!
     PK1_per = 0.0_pReal
     do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
       volfrac = constitutive_matVolFrac(grain,CPFEM_in,cp_en)*constitutive_texVolFrac(grain,CPFEM_in,cp_en)
       PK1_per = PK1_per + volfrac*PK1(:,:,grain)
     enddo
     CPFEM_dPdF_bar(:,:,ip,jp,CPFEM_in,cp_en) = (PK1_per - CPFEM_PK1_bar(:,:,CPFEM_in,cp_en))/1.0e-5_pReal
  enddo
  enddo
 endif
!
 return
!
 END SUBROUTINE
!
!
!********************************************************************
! Calculates the relaxed deformation gradients of grains
!********************************************************************
 subroutine GIA_RelaxedDeformation(&
     F,&        ! relaxed deformation gradient of grains
     F_bar,&    ! effective deformation gradient
     r)         ! relaxation vectors at boundary
!
 implicit none
!
 real(pReal), dimension(3,3)   :: F_bar
 real(pReal), dimension(3,3,8) :: F
 real(pReal), dimension(3,12)  :: r,n
 integer(pInt) i,j,iBoun,grain
!
 n = GIA_bNorm
 do i = 1,3
 do j = 1,3
   F(i,j,1) = F_bar(i,j) + n(i, 1)*r(j, 1) + n(i, 5)*r(j, 5) + n(i, 9)*r(j, 9)
   F(i,j,2) = F_bar(i,j) - n(i, 1)*r(j, 1) + n(i, 6)*r(j, 6) + n(i,10)*r(j,10)
   F(i,j,3) = F_bar(i,j) + n(i, 2)*r(j, 2) - n(i, 5)*r(j, 5) + n(i,11)*r(j,11)
   F(i,j,4) = F_bar(i,j) - n(i, 2)*r(j, 2) - n(i, 6)*r(j, 6) + n(i,12)*r(j,12)
   F(i,j,5) = F_bar(i,j) + n(i, 3)*r(j, 3) + n(i, 7)*r(j, 7) - n(i, 9)*r(j, 9)
   F(i,j,6) = F_bar(i,j) - n(i, 3)*r(j, 3) + n(i, 8)*r(j, 8) - n(i,10)*r(j,10)
   F(i,j,7) = F_bar(i,j) + n(i, 4)*r(j, 4) - n(i, 7)*r(j, 7) - n(i,11)*r(j,11)
   F(i,j,8) = F_bar(i,j) - n(i, 4)*r(j, 4) - n(i, 8)*r(j, 8) - n(i,12)*r(j,12)
 enddo
 enddo
!
 return
!
 END SUBROUTINE
!
!
!********************************************************************
! Calculates the jump of tensors across the grain boundary
!********************************************************************
 subroutine GIA_BoundaryJump(&
     F_boun,&    ! tensor jump across the boundary
     F_bulk)     ! bulk tensor
!
 implicit none
!
 real(pReal), dimension(3,3,12) :: F_boun
 real(pReal), dimension(3,3,8)  :: F_bulk
 integer(pInt) i,j,iBoun,grain
!
 F_boun(:,:, 1) = F_bulk(:,:,2) - F_bulk(:,:,1)
 F_boun(:,:, 2) = F_bulk(:,:,4) - F_bulk(:,:,3)
 F_boun(:,:, 3) = F_bulk(:,:,6) - F_bulk(:,:,5)
 F_boun(:,:, 4) = F_bulk(:,:,8) - F_bulk(:,:,7)
 F_boun(:,:, 5) = F_bulk(:,:,3) - F_bulk(:,:,1)
 F_boun(:,:, 6) = F_bulk(:,:,4) - F_bulk(:,:,2)
 F_boun(:,:, 7) = F_bulk(:,:,7) - F_bulk(:,:,5)
 F_boun(:,:, 8) = F_bulk(:,:,8) - F_bulk(:,:,6)
 F_boun(:,:, 9) = F_bulk(:,:,5) - F_bulk(:,:,1)
 F_boun(:,:,10) = F_bulk(:,:,6) - F_bulk(:,:,2)
 F_boun(:,:,11) = F_bulk(:,:,7) - F_bulk(:,:,3)
 F_boun(:,:,12) = F_bulk(:,:,8) - F_bulk(:,:,4)
!
 return
!
 END SUBROUTINE
!
!
!********************************************************************
! Calculates the jump of tensors across the grain boundary
!********************************************************************
 subroutine GIA_JacobianMatrix(&
     dresdvar,&    ! Jacobian matrix
     dPdF,&        ! stress consistent tangent of bulk
     dRdX)         ! stress-like penalty tangent at boundary
!
 implicit none
!
 real(pReal), dimension(3,3,3,3,8)  :: dPdF
 real(pReal), dimension(3,3,3,3,12) :: dRdX
 real(pReal), dimension(36,36)      :: dresdvar
 real(pReal), dimension(3,12)       :: n
 integer(pInt) i,j,k,l
!
 n = GIA_bNorm
 dresdvar = 0.0_pReal
 do i = 1,3
 do k = 1,3
   do l = 1,3
   do j = 1,3
!
!     at boundary 1, influenced by boundary +5, -6, +9, -10
   dresdvar(( 1-1)*3 + j,( 1-1)*3 + l) = dresdvar(( 1-1)*3 + j,( 1-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 1) + dPdF(i,j,k,l, 2))*n(i, 1)*n(k, 1) &
                                         + (dRdX(i,j,k,l, 1) + dRdX(i,j,k,l, 1))*n(i, 1)*n(k, 1)
   dresdvar(( 1-1)*3 + j,( 5-1)*3 + l) = dresdvar(( 1-1)*3 + j,( 5-1)*3 + l) + dPdF(i,j,k,l, 1)*n(i, 1)*n(k, 5) &
                                                                             + dRdX(i,j,k,l, 1)*n(i, 1)*n(k, 5)
   dresdvar(( 1-1)*3 + j,( 6-1)*3 + l) = dresdvar(( 1-1)*3 + j,( 6-1)*3 + l) - dPdF(i,j,k,l, 2)*n(i, 1)*n(k, 6) &
                                                                             - dRdX(i,j,k,l, 1)*n(i, 1)*n(k, 6)
   dresdvar(( 1-1)*3 + j,( 9-1)*3 + l) = dresdvar(( 1-1)*3 + j,( 9-1)*3 + l) + dPdF(i,j,k,l, 1)*n(i, 1)*n(k, 9) &
                                                                             + dRdX(i,j,k,l, 1)*n(i, 1)*n(k, 9)
   dresdvar(( 1-1)*3 + j,(10-1)*3 + l) = dresdvar(( 1-1)*3 + j,(10-1)*3 + l) - dPdF(i,j,k,l, 2)*n(i, 1)*n(k,10) &
                                                                             - dRdX(i,j,k,l, 1)*n(i, 1)*n(k,10)
!
!     at boundary 2, influenced by boundary -5, +6, +11, -12
   dresdvar(( 2-1)*3 + j,( 2-1)*3 + l) = dresdvar(( 2-1)*3 + j,( 2-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 3) + dPdF(i,j,k,l, 4))*n(i, 2)*n(k, 2) &
                                         + (dRdX(i,j,k,l, 2) + dRdX(i,j,k,l, 2))*n(i, 2)*n(k, 2)
   dresdvar(( 2-1)*3 + j,( 5-1)*3 + l) = dresdvar(( 2-1)*3 + j,( 5-1)*3 + l) - dPdF(i,j,k,l, 3)*n(i, 2)*n(k, 5) &
                                                                             - dRdX(i,j,k,l, 2)*n(i, 2)*n(k, 5)
   dresdvar(( 2-1)*3 + j,( 6-1)*3 + l) = dresdvar(( 2-1)*3 + j,( 6-1)*3 + l) + dPdF(i,j,k,l, 4)*n(i, 2)*n(k, 6) &
                                                                             + dRdX(i,j,k,l, 2)*n(i, 2)*n(k, 6)
   dresdvar(( 2-1)*3 + j,(11-1)*3 + l) = dresdvar(( 2-1)*3 + j,(11-1)*3 + l) + dPdF(i,j,k,l, 3)*n(i, 2)*n(k,11) &
                                                                             + dRdX(i,j,k,l, 2)*n(i, 2)*n(k,11)
   dresdvar(( 2-1)*3 + j,(12-1)*3 + l) = dresdvar(( 2-1)*3 + j,(12-1)*3 + l) - dPdF(i,j,k,l, 4)*n(i, 2)*n(k,12) &
                                                                             - dRdX(i,j,k,l, 2)*n(i, 2)*n(k,12)
!
!     at boundary 3, influenced by boundary +7, -8, -9, +10
   dresdvar(( 3-1)*3 + j,( 3-1)*3 + l) = dresdvar(( 3-1)*3 + j,( 3-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 5) + dPdF(i,j,k,l, 6))*n(i, 3)*n(k, 3) &
                                         + (dRdX(i,j,k,l, 3) + dRdX(i,j,k,l, 3))*n(i, 3)*n(k, 3)
   dresdvar(( 3-1)*3 + j,( 7-1)*3 + l) = dresdvar(( 3-1)*3 + j,( 7-1)*3 + l) + dPdF(i,j,k,l, 5)*n(i, 3)*n(k, 7) &
                                                                             + dRdX(i,j,k,l, 3)*n(i, 3)*n(k, 7)
   dresdvar(( 3-1)*3 + j,( 8-1)*3 + l) = dresdvar(( 3-1)*3 + j,( 8-1)*3 + l) - dPdF(i,j,k,l, 6)*n(i, 3)*n(k, 8) &
                                                                             - dRdX(i,j,k,l, 3)*n(i, 3)*n(k, 8)
   dresdvar(( 3-1)*3 + j,( 9-1)*3 + l) = dresdvar(( 3-1)*3 + j,( 9-1)*3 + l) - dPdF(i,j,k,l, 5)*n(i, 3)*n(k, 9) &
                                                                             - dRdX(i,j,k,l, 3)*n(i, 3)*n(k, 9)
   dresdvar(( 3-1)*3 + j,(10-1)*3 + l) = dresdvar(( 3-1)*3 + j,(10-1)*3 + l) + dPdF(i,j,k,l, 6)*n(i, 3)*n(k,10) &
                                                                             + dRdX(i,j,k,l, 3)*n(i, 3)*n(k,10)
!
!     at boundary 4, influenced by boundary -7, +8, -11, +12
   dresdvar(( 4-1)*3 + j,( 4-1)*3 + l) = dresdvar(( 4-1)*3 + j,( 4-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 7) + dPdF(i,j,k,l, 8))*n(i, 4)*n(k, 4) &
                                         + (dRdX(i,j,k,l, 4) + dRdX(i,j,k,l, 4))*n(i, 4)*n(k, 4)
   dresdvar(( 4-1)*3 + j,( 7-1)*3 + l) = dresdvar(( 4-1)*3 + j,( 7-1)*3 + l) - dPdF(i,j,k,l, 7)*n(i, 4)*n(k, 7) &
                                                                             - dRdX(i,j,k,l, 4)*n(i, 4)*n(k, 7)
   dresdvar(( 4-1)*3 + j,( 8-1)*3 + l) = dresdvar(( 4-1)*3 + j,( 8-1)*3 + l) + dPdF(i,j,k,l, 8)*n(i, 4)*n(k, 8) &
                                                                             + dRdX(i,j,k,l, 4)*n(i, 4)*n(k, 8)
   dresdvar(( 4-1)*3 + j,(11-1)*3 + l) = dresdvar(( 4-1)*3 + j,(11-1)*3 + l) - dPdF(i,j,k,l, 7)*n(i, 4)*n(k,11) &
                                                                             - dRdX(i,j,k,l, 4)*n(i, 4)*n(k,11)
   dresdvar(( 4-1)*3 + j,(12-1)*3 + l) = dresdvar(( 4-1)*3 + j,(12-1)*3 + l) + dPdF(i,j,k,l, 8)*n(i, 4)*n(k,12) &
                                                                             + dRdX(i,j,k,l, 4)*n(i, 4)*n(k,12)
!
!     at boundary 5, influenced by boundary +1, -2, +9, -11
   dresdvar(( 5-1)*3 + j,( 5-1)*3 + l) = dresdvar(( 5-1)*3 + j,( 5-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 1) + dPdF(i,j,k,l, 3))*n(i, 5)*n(k, 5) &
                                         + (dRdX(i,j,k,l, 5) + dRdX(i,j,k,l, 5))*n(i, 5)*n(k, 5)
   dresdvar(( 5-1)*3 + j,( 1-1)*3 + l) = dresdvar(( 5-1)*3 + j,( 1-1)*3 + l) + dPdF(i,j,k,l, 1)*n(i, 5)*n(k, 1) &
                                                                             + dRdX(i,j,k,l, 5)*n(i, 5)*n(k, 1)
   dresdvar(( 5-1)*3 + j,( 2-1)*3 + l) = dresdvar(( 5-1)*3 + j,( 2-1)*3 + l) - dPdF(i,j,k,l, 3)*n(i, 5)*n(k, 2) &
                                                                             - dRdX(i,j,k,l, 5)*n(i, 5)*n(k, 2)
   dresdvar(( 5-1)*3 + j,( 9-1)*3 + l) = dresdvar(( 5-1)*3 + j,( 9-1)*3 + l) + dPdF(i,j,k,l, 1)*n(i, 5)*n(k, 9) &
                                                                             + dRdX(i,j,k,l, 5)*n(i, 5)*n(k, 9)
   dresdvar(( 5-1)*3 + j,(11-1)*3 + l) = dresdvar(( 5-1)*3 + j,(11-1)*3 + l) - dPdF(i,j,k,l, 3)*n(i, 5)*n(k,11) &
                                                                             - dRdX(i,j,k,l, 5)*n(i, 5)*n(k,11)
!
!     at boundary 6, influenced by boundary -1, +2, +10, -12
   dresdvar(( 6-1)*3 + j,( 6-1)*3 + l) = dresdvar(( 6-1)*3 + j,( 6-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 2) + dPdF(i,j,k,l, 4))*n(i, 6)*n(k, 6) &
                                         + (dRdX(i,j,k,l, 6) + dRdX(i,j,k,l, 6))*n(i, 6)*n(k, 6)
   dresdvar(( 6-1)*3 + j,( 1-1)*3 + l) = dresdvar(( 6-1)*3 + j,( 1-1)*3 + l) - dPdF(i,j,k,l, 2)*n(i, 6)*n(k, 1) &
                                                                             - dRdX(i,j,k,l, 6)*n(i, 6)*n(k, 1)
   dresdvar(( 6-1)*3 + j,( 2-1)*3 + l) = dresdvar(( 6-1)*3 + j,( 2-1)*3 + l) + dPdF(i,j,k,l, 4)*n(i, 6)*n(k, 2) &
                                                                             + dRdX(i,j,k,l, 6)*n(i, 6)*n(k, 2)
   dresdvar(( 6-1)*3 + j,(10-1)*3 + l) = dresdvar(( 6-1)*3 + j,(10-1)*3 + l) + dPdF(i,j,k,l, 2)*n(i, 6)*n(k,10) &
                                                                             + dRdX(i,j,k,l, 6)*n(i, 6)*n(k,10)
   dresdvar(( 6-1)*3 + j,(12-1)*3 + l) = dresdvar(( 6-1)*3 + j,(12-1)*3 + l) - dPdF(i,j,k,l, 4)*n(i, 6)*n(k,12) &
                                                                             - dRdX(i,j,k,l, 6)*n(i, 6)*n(k,12)
!
!     at boundary 7, influenced by boundary +3, -4, -9, +11
   dresdvar(( 7-1)*3 + j,( 7-1)*3 + l) = dresdvar(( 7-1)*3 + j,( 7-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 5) + dPdF(i,j,k,l, 7))*n(i, 7)*n(k, 7) &
                                         + (dRdX(i,j,k,l, 7) + dRdX(i,j,k,l, 7))*n(i, 7)*n(k, 7)
   dresdvar(( 7-1)*3 + j,( 3-1)*3 + l) = dresdvar(( 7-1)*3 + j,( 3-1)*3 + l) + dPdF(i,j,k,l, 5)*n(i, 7)*n(k, 3) &
                                                                             + dRdX(i,j,k,l, 7)*n(i, 7)*n(k, 3)
   dresdvar(( 7-1)*3 + j,( 4-1)*3 + l) = dresdvar(( 7-1)*3 + j,( 4-1)*3 + l) - dPdF(i,j,k,l, 7)*n(i, 7)*n(k, 4) &
                                                                             - dRdX(i,j,k,l, 7)*n(i, 7)*n(k, 4)
   dresdvar(( 7-1)*3 + j,( 9-1)*3 + l) = dresdvar(( 7-1)*3 + j,( 9-1)*3 + l) - dPdF(i,j,k,l, 5)*n(i, 7)*n(k, 9) &
                                                                             - dRdX(i,j,k,l, 7)*n(i, 7)*n(k, 9)
   dresdvar(( 7-1)*3 + j,(11-1)*3 + l) = dresdvar(( 7-1)*3 + j,(11-1)*3 + l) + dPdF(i,j,k,l, 7)*n(i, 7)*n(k,11) &
                                                                             + dRdX(i,j,k,l, 7)*n(i, 7)*n(k,11)
!
!     at boundary 8, influenced by boundary -3, +4, -10, +12
   dresdvar(( 8-1)*3 + j,( 8-1)*3 + l) = dresdvar(( 8-1)*3 + j,( 8-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 6) + dPdF(i,j,k,l, 8))*n(i, 8)*n(k, 8) &
                                         + (dRdX(i,j,k,l, 8) + dRdX(i,j,k,l, 8))*n(i, 8)*n(k, 8)
   dresdvar(( 8-1)*3 + j,( 3-1)*3 + l) = dresdvar(( 8-1)*3 + j,( 3-1)*3 + l) - dPdF(i,j,k,l, 6)*n(i, 8)*n(k, 3) &
                                                                             - dRdX(i,j,k,l, 8)*n(i, 8)*n(k, 3)
   dresdvar(( 8-1)*3 + j,( 4-1)*3 + l) = dresdvar(( 8-1)*3 + j,( 4-1)*3 + l) + dPdF(i,j,k,l, 8)*n(i, 8)*n(k, 4) &
                                                                             + dRdX(i,j,k,l, 8)*n(i, 8)*n(k, 4)
   dresdvar(( 8-1)*3 + j,(10-1)*3 + l) = dresdvar(( 8-1)*3 + j,(10-1)*3 + l) - dPdF(i,j,k,l, 6)*n(i, 8)*n(k,10) &
                                                                             - dRdX(i,j,k,l, 8)*n(i, 8)*n(k,10)
   dresdvar(( 8-1)*3 + j,(12-1)*3 + l) = dresdvar(( 8-1)*3 + j,(12-1)*3 + l) + dPdF(i,j,k,l, 8)*n(i, 8)*n(k,12) &
                                                                             + dRdX(i,j,k,l, 8)*n(i, 8)*n(k,12)
!
!     at boundary 9, influenced by boundary +1, -3, +5, -7
   dresdvar(( 9-1)*3 + j,( 9-1)*3 + l) = dresdvar(( 9-1)*3 + j,( 9-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 1) + dPdF(i,j,k,l, 5))*n(i, 9)*n(k, 9) &
                                         + (dRdX(i,j,k,l, 9) + dRdX(i,j,k,l, 9))*n(i, 9)*n(k, 9)
   dresdvar(( 9-1)*3 + j,( 1-1)*3 + l) = dresdvar(( 9-1)*3 + j,( 1-1)*3 + l) + dPdF(i,j,k,l, 1)*n(i, 9)*n(k, 1) &
                                                                             + dRdX(i,j,k,l, 9)*n(i, 9)*n(k, 1)
   dresdvar(( 9-1)*3 + j,( 3-1)*3 + l) = dresdvar(( 9-1)*3 + j,( 3-1)*3 + l) - dPdF(i,j,k,l, 5)*n(i, 9)*n(k, 3) &
                                                                             - dRdX(i,j,k,l, 9)*n(i, 9)*n(k, 3)
   dresdvar(( 9-1)*3 + j,( 5-1)*3 + l) = dresdvar(( 9-1)*3 + j,( 5-1)*3 + l) + dPdF(i,j,k,l, 1)*n(i, 9)*n(k, 5) &
                                                                             + dRdX(i,j,k,l, 9)*n(i, 9)*n(k, 5)
   dresdvar(( 9-1)*3 + j,( 7-1)*3 + l) = dresdvar(( 9-1)*3 + j,( 7-1)*3 + l) - dPdF(i,j,k,l, 5)*n(i, 9)*n(k, 7) &
                                                                             - dRdX(i,j,k,l, 9)*n(i, 9)*n(k, 7)
!
!     at boundary 10, influenced by boundary -1, +3, +6, -8
   dresdvar((10-1)*3 + j,(10-1)*3 + l) = dresdvar((10-1)*3 + j,(10-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 2) + dPdF(i,j,k,l, 6))*n(i,10)*n(k,10) &
                                         + (dRdX(i,j,k,l,10) + dRdX(i,j,k,l,10))*n(i,10)*n(k,10)
   dresdvar((10-1)*3 + j,( 1-1)*3 + l) = dresdvar((10-1)*3 + j,( 1-1)*3 + l) - dPdF(i,j,k,l, 2)*n(i,10)*n(k, 1) &
                                                                             - dRdX(i,j,k,l,10)*n(i,10)*n(k, 1)
   dresdvar((10-1)*3 + j,( 3-1)*3 + l) = dresdvar((10-1)*3 + j,( 3-1)*3 + l) + dPdF(i,j,k,l, 6)*n(i,10)*n(k, 3) &
                                                                             + dRdX(i,j,k,l,10)*n(i,10)*n(k, 3)
   dresdvar((10-1)*3 + j,( 6-1)*3 + l) = dresdvar((10-1)*3 + j,( 6-1)*3 + l) + dPdF(i,j,k,l, 2)*n(i,10)*n(k, 6) &
                                                                             + dRdX(i,j,k,l,10)*n(i,10)*n(k, 6)
   dresdvar((10-1)*3 + j,( 8-1)*3 + l) = dresdvar((10-1)*3 + j,( 8-1)*3 + l) - dPdF(i,j,k,l, 6)*n(i,10)*n(k, 8) &
                                                                             - dRdX(i,j,k,l,10)*n(i,10)*n(k, 8)
!
!     at boundary 11, influenced by boundary +2, -4, -5, +7
   dresdvar((11-1)*3 + j,(11-1)*3 + l) = dresdvar((11-1)*3 + j,(11-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 3) + dPdF(i,j,k,l, 7))*n(i,11)*n(k,11) &
                                         + (dRdX(i,j,k,l,11) + dRdX(i,j,k,l,11))*n(i,11)*n(k,11)
   dresdvar((11-1)*3 + j,( 2-1)*3 + l) = dresdvar((11-1)*3 + j,( 2-1)*3 + l) + dPdF(i,j,k,l, 3)*n(i,11)*n(k, 2) &
                                                                             + dRdX(i,j,k,l,11)*n(i,11)*n(k, 2)
   dresdvar((11-1)*3 + j,( 4-1)*3 + l) = dresdvar((11-1)*3 + j,( 4-1)*3 + l) - dPdF(i,j,k,l, 7)*n(i,11)*n(k, 4) &
                                                                             - dRdX(i,j,k,l,11)*n(i,11)*n(k, 4)
   dresdvar((11-1)*3 + j,( 5-1)*3 + l) = dresdvar((11-1)*3 + j,( 5-1)*3 + l) - dPdF(i,j,k,l, 3)*n(i,11)*n(k, 5) &
                                                                             - dRdX(i,j,k,l,11)*n(i,11)*n(k, 5)
   dresdvar((11-1)*3 + j,( 7-1)*3 + l) = dresdvar((11-1)*3 + j,( 7-1)*3 + l) + dPdF(i,j,k,l, 7)*n(i,11)*n(k, 7) &
                                                                             + dRdX(i,j,k,l,11)*n(i,11)*n(k, 7)
!
!     at boundary 12, influenced by boundary -2, +4, -6, +8
   dresdvar((12-1)*3 + j,(12-1)*3 + l) = dresdvar((12-1)*3 + j,(12-1)*3 + l) &
                                         + (dPdF(i,j,k,l, 4) + dPdF(i,j,k,l, 8))*n(i,12)*n(k,12) &
                                         + (dRdX(i,j,k,l,12) + dRdX(i,j,k,l,12))*n(i,12)*n(k,12)
   dresdvar((12-1)*3 + j,( 2-1)*3 + l) = dresdvar((12-1)*3 + j,( 2-1)*3 + l) - dPdF(i,j,k,l, 4)*n(i,12)*n(k, 2) &
                                                                             - dRdX(i,j,k,l,12)*n(i,12)*n(k, 2)
   dresdvar((12-1)*3 + j,( 4-1)*3 + l) = dresdvar((12-1)*3 + j,( 4-1)*3 + l) + dPdF(i,j,k,l, 8)*n(i,12)*n(k, 4) &
                                                                             + dRdX(i,j,k,l,12)*n(i,12)*n(k, 4)
   dresdvar((12-1)*3 + j,( 6-1)*3 + l) = dresdvar((12-1)*3 + j,( 6-1)*3 + l) - dPdF(i,j,k,l, 4)*n(i,12)*n(k, 6) &
                                                                             - dRdX(i,j,k,l,12)*n(i,12)*n(k, 6)
   dresdvar((12-1)*3 + j,( 8-1)*3 + l) = dresdvar((12-1)*3 + j,( 8-1)*3 + l) + dPdF(i,j,k,l, 8)*n(i,12)*n(k, 8) &
                                                                             + dRdX(i,j,k,l,12)*n(i,12)*n(k, 8)
!
   enddo
   enddo
 enddo
 enddo
!
 return
!
 END SUBROUTINE
!
!
 END MODULE
!##############################################################

