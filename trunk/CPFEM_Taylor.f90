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
 real(pReal), dimension (:,:,:,:,:),  allocatable :: CPFEM_Lp
 real(pReal), dimension (:,:,:,:,:),  allocatable :: CPFEM_Fp_old
 real(pReal), dimension (:,:,:,:,:),  allocatable :: CPFEM_Fp_new
 real(pReal), parameter :: CPFEM_odd_stress = 1e15_pReal, CPFEM_odd_jacobian = 1e50_pReal
 integer(pInt) :: CPFEM_Nresults   = 4_pInt    ! three Euler angles plus volume fraction
 logical :: CPFEM_init_done = .false.          ! remember if init has been done already
 logical :: CPFEM_calc_done = .false.          ! remember if first IP has already calced the results
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
 integer(pInt) e,i,g
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
!    *** Plastic velocity gradient ***
 allocate(CPFEM_Lp(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Lp = 0.0_pReal

!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fp_new = 0.0_pReal
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall (e=1:mesh_NcpElems,i=1:mesh_maxNips,g=1:constitutive_maxNgrains) &
   CPFEM_Fp_old(:,:,g,i,e) = math_EulerToR(constitutive_EulerAngles(:,g,i,e))  ! plastic def gradient reflects init orientation
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
 write(6,*) 'CPFEM_Lp:            ', shape(CPFEM_Lp)
 write(6,*) 'CPFEM_Fp_old:        ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:        ', shape(CPFEM_Fp_new)
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
 use math, only: math_init, invnrmMandel, math_identity2nd, math_Mandel3333to66,math_Mandel33to6,math_Mandel6to33,math_det3x3,math_I3
 use mesh, only: mesh_init,mesh_FEasCP, mesh_NcpElems, FE_Nips, FE_mapElemtype, mesh_element
 use lattice, only: lattice_init
 use constitutive, only: constitutive_init,constitutive_state_old,constitutive_state_new,material_Cslip_66
 implicit none
!
 integer(pInt) CPFEM_en, CPFEM_in, cp_en, CPFEM_ngens, i,j,k,l,m,n, e
 real(pReal), dimension (3,3)        :: ffn,ffn1,Kirchhoff_bar
 real(pReal), dimension (3,3,3,3)    :: H_bar
 real(pReal), dimension(CPFEM_ngens) :: CPFEM_stress
 real(pReal), dimension(CPFEM_ngens,CPFEM_ngens) :: CPFEM_jaco, odd_jaco
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
  if (cp_en == 1 .and. CPFEM_in == 1) &
    write(6,'(a10,1x,f8.4,1x,a10,1x,i4,1x,a10,1x,i3,1x,a10,1x,i2,x,a10,1x,i2)') &
    'theTime',theTime,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,&
    'mode',CPFEM_mode
!
 select case (CPFEM_mode)
    case (2,1)     ! regular computation (with aging of results)
       if (any(abs(ffn1 - CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en)) > relevantStrain)) then
           CPFEM_stress_bar(1:CPFEM_ngens,:,:) = CPFEM_odd_stress
           odd_jaco = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
           forall (i=1:mesh_NcpElems, j=1:FE_Nips(mesh_element(2,e))) &
              CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,j,i) = odd_jaco
           outdatedFFN1 = .true.
           if (.not. CPFEM_calc_done .AND.CPFEM_mode == 1) then
             CPFEM_Fp_old            = CPFEM_Fp_new
             constitutive_state_old  = constitutive_state_new
           endif
!$OMP CRITICAL (write2out)
           write(6,*) 'WARNING: FFN1 changed for ip', CPFEM_in, 'of element', cp_en
!$OMP END CRITICAL (write2out)
           return
       else
         if (.not. CPFEM_calc_done) then                ! puuh, me needs doing all the work...
           if (CPFEM_mode == 1) then                  ! age results at start of new increment
             CPFEM_Fp_old            = CPFEM_Fp_new
             constitutive_state_old  = constitutive_state_new
           endif
           debug_cutbackDistribution   = 0_pInt       ! initialize debugging data
           debug_InnerLoopDistribution = 0_pInt
           debug_OuterLoopDistribution = 0_pInt
!
!$OMP PARALLEL DO
           do e=1,mesh_NcpElems                       ! ## this shall be done in a parallel loop in the future ##
               do i=1,FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
                  debugger = (e==1 .and. i==1)        ! switch on debugging for first IP in first element
                  call CPFEM_MaterialPoint(CPFEM_updateJaco, CPFEM_dt, i, e)
               enddo
           enddo
!$OMP END PARALLEL DO
           call debug_info()                          ! output of debugging/performance statistics
           CPFEM_calc_done = .true.                   ! now calc is done
         endif    
!       translate from P and dP/dF to CS and dCS/dE
!$OMP CRITICAL (evilmatmul)
         Kirchhoff_bar = matmul(CPFEM_PK1_bar(:,:,CPFEM_in, cp_en),transpose(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en)))
!$OMP END CRITICAL (evilmatmul)
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
         CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel3333to66(J_inverse*H_bar)
! if (CPFEM_in==8 .and. cp_en==80) then
!   do e=1,80
!       do i=1,8
!       write(6,*)
!       write(6,*) e, i
!       write(6,*) CPFEM_stress_bar(1:CPFEM_ngens,i,e)
!       write(6,*)
!       write(6,*) CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,i,e)
!       enddo
!   enddo
! endif
!
      endif
    case (3)    ! collect and return odd result
       CPFEM_Temperature(CPFEM_in,cp_en)  = Temperature
       CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)  = ffn
       CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) = ffn1
       CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_stress
       CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       CPFEM_calc_done = .false.
! if (CPFEM_in==8 .and. cp_en==80) then
!   do e=1,80
!       do i=1,8
!       write(6,*)
!       write(6,*) e, i
!       write(6,*) ffn1
!       enddo
!   enddo
! endif

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
 integer(pInt) cp_en,CPFEM_in,grain
 logical updateJaco,error
 real(pReal) CPFEM_dt,volfrac
 real(pReal), dimension(3,3)        :: U,R,Fe1
 real(pReal), dimension(3,3)        :: PK1
 real(pReal), dimension(3,3,3,3)    :: dPdF,dPdF_bar_old
!
 CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) = 0.0_pReal                         ! zero out average first PK stress
 if (updateJaco) then
   dPdF_bar_old = CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en)               ! remember former average consistent tangent
   CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) = 0.0_pReal                  ! zero out avg consistent tangent for later assembly
 endif

 do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
   dPdF = dPdF_bar_old                                                 ! preguess consistent tangent of grain with avg
   call SingleCrystallite(msg,PK1,dPdF,&
                      CPFEM_results(5:4+constitutive_Nresults(grain,CPFEM_in,cp_en),grain,CPFEM_in,cp_en),&
                      CPFEM_Lp(:,:,grain,CPFEM_in,cp_en),&
                      CPFEM_Fp_new(:,:,grain,CPFEM_in,cp_en),Fe1,constitutive_state_new(:,grain,CPFEM_in,cp_en),&   ! output up to here
                      CPFEM_dt,cp_en,CPFEM_in,grain,updateJaco,&
                      CPFEM_Temperature(CPFEM_in,cp_en),&
                      CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en),CPFEM_ffn_bar(:,:,CPFEM_in,cp_en),&
                      CPFEM_Fp_old(:,:,grain,CPFEM_in,cp_en),constitutive_state_old(:,grain,CPFEM_in,cp_en))

   if (msg /= 'ok') then                                               ! solution not reached --> exit
!$OMP CRITICAL (write2out)
       write(6,*) 'grain loop failed to converge @ EL:',cp_en,' IP:',CPFEM_in
!$OMP END CRITICAL (write2out)
       call IO_error(600)
       return 
   endif

   volfrac = constitutive_matVolFrac(grain,CPFEM_in,cp_en)*constitutive_texVolFrac(grain,CPFEM_in,cp_en)
   CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) = CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) + volfrac*PK1
   if (updateJaco) CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) = &
                   CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) + volfrac*dPdF   ! add up crystallite stiffnesses
				                                                           ! (may have "holes" corresponding 
				                                                           ! to former avg tangent)
!
!   update results plotted in MENTAT
   call math_pDecomposition(Fe1,U,R,error) ! polar decomposition
   if (error) then
!$OMP CRITICAL (write2out)
     write(6,*) 'polar decomposition of', Fe1
     write(6,*) 'Grain:             ',grain
     write(6,*) 'Integration point: ',CPFEM_in
     write(6,*) 'Element:           ',mesh_element(1,cp_en)
!$OMP END CRITICAL (write2out)
     call IO_error(650)
     return
   endif
   CPFEM_results(1:3,grain,CPFEM_in,cp_en) = math_RtoEuler(transpose(R))*inDeg        ! orientation
   CPFEM_results(4  ,grain,CPFEM_in,cp_en) = volfrac                                  ! volume fraction of orientation
 enddo ! grain
!
 return
!
 END SUBROUTINE
!
 END MODULE
!##############################################################

