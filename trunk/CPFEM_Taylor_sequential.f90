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
 real(pReal), dimension (:,:,:,:,:),    allocatable :: CPFEM_Fe1
 real(pReal), dimension (:,:,:,:),      allocatable :: CPFEM_Tstar_v
 real(pReal), parameter :: CPFEM_odd_stress = 1e15_pReal, CPFEM_odd_jacobian = 1e50_pReal
 integer(pInt) :: CPFEM_Nresults   = 4_pInt    ! three Euler angles plus volume fraction
 logical :: CPFEM_init_done    = .false.       ! remember whether init has been done already
 logical :: CPFEM_calc_done    = .false.       ! remember whether first IP has already calced the results
 logical :: CPFEM_results_aged = .false.       ! remember whether results have been aged at inc start
!    *** Solution at single crystallite level ***
!
 logical, dimension (:,:,:),allocatable :: crystallite_converged            !individual covergence flag per grain
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
 allocate(CPFEM_Temperature(mesh_maxNips,mesh_NcpElems)) ; CPFEM_Temperature = Temperature
 allocate(CPFEM_ffn_bar(3,3,mesh_maxNips,mesh_NcpElems))
 forall(e=1:mesh_NcpElems,i=1:mesh_maxNips) CPFEM_ffn_bar(:,:,i,e)             = math_I3
 allocate(CPFEM_ffn(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall(g=1:constitutive_maxNgrains,e=1:mesh_NcpElems,i=1:mesh_maxNips) CPFEM_ffn(:,:,g,i,e) = math_I3
 allocate(CPFEM_ffn1_bar(3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1_bar   = CPFEM_ffn_bar
 allocate(CPFEM_ffn1(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1 = CPFEM_ffn
 allocate(CPFEM_PK1_bar(3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_PK1_bar    = 0.0_pReal
 allocate(CPFEM_PK1(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_PK1 = 0.0_pReal
 allocate(CPFEM_dPdF_bar(3,3,3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF_bar = 0.0_pReal
 allocate(CPFEM_dPdF_bar_old(3,3,3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF_bar_old = 0.0_pReal
 allocate(CPFEM_dPdF(3,3,3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_dPdF = 0.0_pReal
 allocate(CPFEM_stress_bar(6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_stress_bar = 0.0_pReal
 allocate(CPFEM_jaco_bar(6,6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_jaco_bar   = 0.0_pReal
 allocate(CPFEM_jaco_knownGood(6,6,mesh_maxNips,mesh_NcpElems)) ;   CPFEM_jaco_knownGood = 0.0_pReal
!
!    *** User defined results !!! MISSING incorporate consti_Nresults ***
 allocate(CPFEM_results(CPFEM_Nresults+constitutive_maxNresults,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_results = 0.0_pReal
!
!    *** Plastic velocity gradient ***
 allocate(CPFEM_Lp_old(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Lp_old = 0.0_pReal
 allocate(CPFEM_Lp_new(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Lp_new = 0.0_pReal

!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fp_new = 0.0_pReal
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall (e=1:mesh_NcpElems,i=1:mesh_maxNips,g=1:constitutive_maxNgrains) &
   CPFEM_Fp_old(:,:,g,i,e) = math_EulerToR(constitutive_EulerAngles(:,g,i,e))  ! plastic def gradient reflects init orientation
!    *** Elastic deformation gradient at (t=t1) ***
 allocate(CPFEM_Fe1(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fe1 = 0.0_pReal
!    *** Stress vector at (t=t1) ***
 allocate(CPFEM_Tstar_v(6,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Tstar_v = 0.0_pReal
!
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
 write(6,*) 'CPFEM_Fe1:            ', shape(CPFEM_Fe1)
 write(6,*) 'CPFEM_Tstar_v:        ', shape(CPFEM_Tstar_v)
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
   call crystallite_init()
   call CPFEM_init(Temperature)			
   CPFEM_init_done = .true.
 endif
!
 cp_en = mesh_FEasCP('elem',CPFEM_en)
 if (cp_en == 1 .and. CPFEM_in == 1) &
   write(6,'(a10,1x,f8.4,1x,a10,1x,i4,1x,a10,1x,i3,1x,a10,1x,i2,x,a10,1x,i2)') &
   'theTime',theTime,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,&
   'mode',CPFEM_mode

 if (CPFEM_mode /= 1) CPFEM_results_aged = .false.
 
 select case (CPFEM_mode)

   case (2,1)                                     ! *** regular computation (with aging of results) ***
     if (CPFEM_mode == 1 .and. &
		 .not. CPFEM_results_aged) then           ! age results at start of new increment
         CPFEM_Fp_old            = CPFEM_Fp_new
         constitutive_state_old  = constitutive_state_new
         CPFEM_results_aged = .true.              ! aging is done
		 write (6,*) ')))))))))))))) results aged (((((((((((((((',cp_en,CPFEM_in
     endif

     CPFEM_Temperature(CPFEM_in,cp_en)  = Temperature   ! store temperature
     CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)  = ffn           ! store def grad for start of inc
     CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) = ffn1          ! store def grad for end of inc
     debugger = (cp_en == 1160 .and. CPFEM_in == 6)      ! switch on debugging
     call CPFEM_MaterialPoint(CPFEM_updateJaco, CPFEM_dt, CPFEM_in, cp_en)   ! call for result at this IP

! translate from P and dP/dF to CS and dCS/dE

     Kirchhoff_bar = math_mul33x33(CPFEM_PK1_bar(:,:,CPFEM_in, cp_en),transpose(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en)))
     J_inverse  = 1.0_pReal/math_det3x3(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en))
     CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff_bar)

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

   case (3)    ! *** collect and return odd result ***
     CPFEM_Temperature(CPFEM_in,cp_en)  = Temperature
     CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)  = ffn
     CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) = ffn1
     CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_stress
     CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
     CPFEM_calc_done = .false.

   case (4)    ! *** do nothing since we can recycle the former results (MARC specialty) ***
   case (5)    ! *** record consistent tangent at beginning of new increment ***
     CPFEM_jaco_knownGood = CPFEM_jaco_bar
   case (6)    ! *** restore consistent tangent after cutback ***
     CPFEM_jaco_bar = CPFEM_jaco_knownGood
 end select
!
! return the local stress and the jacobian from storage
 CPFEM_stress(1:CPFEM_ngens) = CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en)
 CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en)
 if (debugger) write (6,'(a,/,6(6(f9.3,x)/))') 'stiffness / GPa',CPFEM_jaco(1:CPFEM_ngens,:)/1e9_pReal
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
 use mesh, only: mesh_element, mesh_NcpElems, FE_Nips
! use crystallite
 use constitutive
 implicit none

!
 integer(pInt) cp_en,CPFEM_in,g,i,e
 integer(pInt) el_start, el_end, ip_start, ip_end
 logical updateJaco,error
 real(pReal) CPFEM_dt,volfrac
 real(pReal), dimension(3,3)        :: U,R !,Fe1
! real(pReal), dimension(3,3)        :: PK1
! real(pReal), dimension(3,3,3,3)    :: dPdF,dPdF_bar_old
!
 CPFEM_PK1_bar = 0.0_pReal                         ! zero out average first PK stress
!initialize element loop
 if (cp_en /= 0_pInt) then
    el_start = cp_en
    el_end = cp_en
 else
    el_start = 1_pInt
    el_end = mesh_NcpElems
 endif
! prescribe FFN and FFN1 depending on homogenization scheme
!$OMP PARALLEL DO
 do e=el_start,el_end
    if(CPFEM_in /= 0_pInt) then
        ip_start = CPFEM_in
        ip_end = CPFEM_in
    else
        ip_start = 1
        ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
    endif
    do i=ip_start,ip_end
        do g=1,texture_Ngrains(mesh_element(4,e))
            CPFEM_ffn(:,:,g,i,e) = CPFEM_ffn_bar(:,:,i,e)    !Taylor homogenization
            CPFEM_ffn1(:,:,g,i,e) = CPFEM_ffn1_bar(:,:,i,e)  !Taylor homogenization
        end do
    end do
 end do
!$OMP END PARALLEL DO
! calculate stress, update state and update jacobian in case needed for all or one ip
 if (updateJaco) then
   CPFEM_dPdF_bar_old = CPFEM_dPdF_bar               ! remember former average consistent tangent
   CPFEM_dPdF_bar = 0.0_pReal                        ! zero out avg consistent tangent for later assembly
 endif
 call SingleCrystallite(updateJaco,CPFEM_dt,el_start,el_end,CPFEM_in)
!******************************************************************************************************
! check convergence of homogenization in case needed
!******************************************************************************************************
! calculate average quantities per ip and post results
!$OMP PARALLEL DO
 do e=el_start,el_end
    if(CPFEM_in /= 0_pInt) then
        ip_start = CPFEM_in
        ip_end = CPFEM_in
    else
        ip_start = 1
        ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
    endif
    do i=ip_start,ip_end
        do g=1,texture_Ngrains(mesh_element(4,e))
            volfrac = constitutive_matVolFrac(g,i,e)*constitutive_texVolFrac(g,i,e)
            CPFEM_PK1_bar(:,:,i,e) = CPFEM_PK1_bar(:,:,i,e) + volfrac * CPFEM_PK1(:,:,g,i,e)
            if (updateJaco) CPFEM_dPdF_bar(:,:,:,:,i,e) = &
                   CPFEM_dPdF_bar(:,:,:,:,i,e) + volfrac * CPFEM_dPdF(:,:,:,:,g,i,e)   ! add up crystallite stiffnesses
				                                                                       ! (may have "holes" corresponding 
				                                                                       ! to former avg tangent)
!   update results plotted in MENTAT
            call math_pDecomposition(CPFEM_Fe1(:,:,g,i,e),U,R,error) ! polar decomposition
            if (error) then
!$OMP CRITICAL (write2out)
                write(6,*) 'polar decomposition of', CPFEM_Fe1(:,:,g,i,e)
                write(6,*) 'Grain:             ',g
                write(6,*) 'Integration point: ',i
                write(6,*) 'Element:           ',mesh_element(1,e)
!$OMP END CRITICAL (write2out)
                call IO_error(650)
                return
            endif
            CPFEM_results(1:3,g,i,e) = math_RtoEuler(transpose(R))*inDeg        ! orientation
            CPFEM_results(4  ,g,i,e) = volfrac                                  ! volume fraction of orientation
        end do
    end do
 end do
!$OMP END PARALLEL DO
!
 return
!
 END SUBROUTINE
!
!
!********************************************************************
! Initialize crystallite
!********************************************************************
 subroutine crystallite_init()
 use mesh, only: mesh_maxNips,mesh_NcpElems
 use constitutive, only: constitutive_maxNgrains

 implicit none

 allocate(crystallite_converged(constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)); crystallite_converged = .false.
!
!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) 'crystallite Initialization'
 write(6,*)
 write(6,*) 'crystallite_converged:      ', shape(crystallite_converged)
 write(6,*)
 call flush(6)
!$OMP END CRITICAL (write2out)
 return
!
 end subroutine
!
!
!********************************************************************
! Calculates the stress and jacobi (if wanted) for all or a single component
!********************************************************************
 subroutine SingleCrystallite(&
     updateJaco,& ! update of Jacobian required
     dt,&         ! time increment
     el_start,&   ! first element in element loop
     el_end,&     ! last element in element loop
     CPFEM_in)    ! IP number
!
 use prec, only: pReal,pInt,pert_Fg,subStepMin, nCutback
 use debug
 use constitutive
 use mesh, only: mesh_element, FE_Nips
 use math
 use IO,   only: IO_error
! use CPFEM

 implicit none
!
 logical updateJaco, JacoOK
 real(preal) dt
 real(pReal), dimension(3,3) :: Fg_pert,Lp_pert, P_pert, Fp_pert, Fe_pert
 real(pReal), dimension(6) :: Tstar_v     
 real(pReal), dimension(constitutive_maxNstatevars) :: state_pert
 integer(pInt) el_start, el_end, CPFEM_in, ip_start, ip_end, g, i, e, k, l, iOuter
!
 crystallite_converged=.true.
!$OMP PARALLEL DO
    do e=el_start,el_end
        if(CPFEM_in /= 0_pInt) then
            ip_start = CPFEM_in
            ip_end = CPFEM_in
        else
            ip_start = 1
            ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
        endif
        do i=ip_start,ip_end
            do g=1,texture_Ngrains(mesh_element(4,e))
                crystallite_converged(g,i,e)=.false.
            end do
        end do
    end do
!$OMP END PARALLEL DO
 constitutive_state_new=constitutive_state_old
 CPFEM_Lp_new = CPFEM_Lp_old
 iOuter = 0_pInt
 do while(any(crystallite_converged(:,:,el_start:el_end))==.false.)
!$OMP PARALLEL DO
    do e=el_start,el_end
        if(CPFEM_in /= 0_pInt) then
            ip_start = CPFEM_in
            ip_end = CPFEM_in
        else
            ip_start = 1
            ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
        endif
        do i=ip_start,ip_end
            do g=1,texture_Ngrains(mesh_element(4,e))
                if(.not.crystallite_converged(g,i,e))&
                    call IntegrateStress(CPFEM_Tstar_v(:,g,i,e), CPFEM_PK1(:,:,g,i,e), CPFEM_ffn1(:,:,g,i,e),&
                     CPFEM_Fp_new(:,:,g,i,e), CPFEM_Fe1(:,:,g,i,e), CPFEM_Lp_new(:,:,g,i,e),&
                     constitutive_state_new(:,g,i,e), dt, g, i, e)
            end do
        end do
    end do
!$OMP END PARALLEL DO
!$OMP PARALLEL DO
    do e=el_start,el_end
        if(CPFEM_in /= 0_pInt) then
            ip_start = CPFEM_in
            ip_end = CPFEM_in
        else
            ip_start = 1
            ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
        endif
        do i=ip_start,ip_end
            do g=1,texture_Ngrains(mesh_element(4,e))
                if(.not.crystallite_converged(g,i,e))& 
                    call UpdateState(CPFEM_Tstar_v(:,g,i,e),constitutive_state_new(:,g,i,e),dt,g,i,e)
            end do
        end do
    end do
!$OMP END PARALLEL DO
    iOuter = iOuter + 1_pInt
    if (iOuter==Nouter) then
!$OMP CRITICAL (write2out)
        write (6,*) 'Terminated outer loop at el,ip,grain',e,i,g
!$OMP CRITICAL (out)
        debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!$OMP END CRITICAL (out)
        call IO_error(600)
!$OMP END CRITICAL (write2out)
    endif
 end do
!$OMP CRITICAL (out)
 debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!$OMP END CRITICAL (out)

 if (wantsConstitutiveResults) then                                ! get the post_results upon request
!$OMP PARALLEL DO
    do e=el_start,el_end
        if(CPFEM_in /= 0_pInt) then
            ip_start = CPFEM_in
            ip_end = CPFEM_in
        else
            ip_start = 1
            ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
        endif
        do i=ip_start,ip_end
            do g=1,texture_Ngrains(mesh_element(4,e))
                CPFEM_results(CPFEM_Nresults+1:CPFEM_Nresults+constitutive_Nresults(g,i,e),g,i,e) =&
                 constitutive_post_results(CPFEM_Tstar_v(:,g,i,e),constitutive_state_new(:,g,i,e),&
                 CPFEM_Temperature(i,e),dt,g,i,e)
            end do
        end do
    end do
!$OMP END PARALLEL DO
 endif
!
!***** Calculate Jacobian *****
 if(updateJaco) then
    if (debugger) then
!$OMP CRITICAL (write2out)
        write (6,*) 'Jacobian calc'
!$OMP END CRITICAL (write2out)
    endif   
!    crystallite_converged=.false.
!$OMP PARALLEL DO
    do e=el_start,el_end
        if(CPFEM_in /= 0_pInt) then
            ip_start = CPFEM_in
            ip_end = CPFEM_in
        else
            ip_start = 1
            ip_end = FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
        endif
        do i=ip_start,ip_end
            do g=1,texture_Ngrains(mesh_element(4,e))
                do k=1,3
                    do l=1,3
                        crystallite_converged(g,i,e)=.false.
                        JacoOK=.true.
                        Fg_pert = CPFEM_ffn1(:,:,g,i,e)                                   ! initialize perturbed Fg
                        Fg_pert(k,l) = Fg_pert(k,l) + pert_Fg                             ! perturb single component
                        Lp_pert = CPFEM_Lp_new(:,:,g,i,e)                                 ! initialize Lp
                        Fp_pert = CPFEM_Fp_new(:,:,g,i,e)                                 ! initialize Fp
                        state_pert = constitutive_state_new(:,g,i,e)                      ! initial guess from end of time step
                        iOuter=0_pInt
                        do while(.not.crystallite_converged(g,i,e))
                            call IntegrateStress(Tstar_v, P_pert, Fg_pert, Fp_pert, Fe_pert, Lp_pert, state_pert, dt, g, i, e)
                            call UpdateState(Tstar_v,state_pert,dt,g,i,e)
                            iOuter = iOuter + 1_pInt
                            if (iOuter==Nouter) then
                                JacoOK=.false.
                                exit
                            endif
                        end do
!$OMP CRITICAL (out)
                        debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!$OMP END CRITICAL (out)
                        if (JacoOK) &
                        CPFEM_dPdF(:,:,k,l,g,i,e) = (P_pert-CPFEM_PK1(:,:,g,i,e))/pert_Fg ! constructing tangent dP_ij/dFg_kl only if valid forward difference
                                                                                          ! otherwise leave component unchanged
                     end do
                end do
            end do
        end do
    end do
!$OMP END PARALLEL DO
 endif
!
 return
!
 end subroutine
!
!********************************************************************
! Update the state for a single component
!********************************************************************
 subroutine UpdateState(&
   Tstar_v,&        ! stress
   state,&          ! state  
   dt,&             ! time increment
   g,&              ! grain number
   i,&              ! integration point number
   e&               ! element number
 )
 use prec, only: pReal,pInt,reltol_Outer
 use constitutive, only: constitutive_dotState, constitutive_state_old, constitutive_Nstatevars
! use CPFEM, only: CPFEM_Temperature
!
 integer(pInt) g, i, e
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_Nstatevars(g, i, e)) :: state, ROuter
 real(pReal) dt
!
 ROuter = state - constitutive_state_old(:,g,i,e) - &
                dt*constitutive_dotState(Tstar_v,state,CPFEM_Temperature(i,e),&
                                         g,i,e)          ! residuum from evolution of microstructure
 state = state - ROuter                                           ! update of microstructure
 if (maxval(abs(ROuter/state),state /= 0.0_pReal) < reltol_Outer) crystallite_converged(g,i,e) = .true.
!
 return
!
 end subroutine
!
!
!********************************************************************
! Calculates the stress for a single component
!********************************************************************
!***********************************************************************
!***     calculation of stress (P), stiffness (dPdF),                ***
!***     and announcment of any                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 subroutine IntegrateStress(&
     Tstar_v,&    ! Stress vector
     P,&          ! first PK stress
     Fg_new,&     ! new global deformation gradient
     Fp_new,&     ! new plastic deformation gradient
     Fe_new,&     ! new "elastic" deformation gradient
     Lp,&         ! plastic velocity gradient
     state_new,&  ! new state variable array
     dt,&         ! time increment
     g,&          ! grain number
     i,&          ! integration point number
     e)           ! element number

!     post_results,& ! plot results from constitutive model
!     Fp_new,&     ! new plastic deformation gradient
!     updateJaco,& ! update of Jacobian required
!     Temperature,& ! temperature of crystallite
!     Fg_old,&     ! old global deformation gradient
!     Fp_old,&     ! old plastic deformation gradient
!     state_old)   ! old state variable array
!
 use prec, only: pReal,pInt,pert_Fg,subStepMin, nCutback
 use debug
 use constitutive, only: constitutive_Nstatevars,constitutive_Nresults,constitutive_state_old
 use math
! use CPFEM
!
 implicit none
!
 character(len=128) msg
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
 real(pReal), dimension(constitutive_Nstatevars(g,i,e)) :: state_new
! real(pReal), dimension(constitutive_Nstatevars(g,i,e)) :: state_current
! 
! debugger= e==1.and.i==1
 deltaFg = Fg_new - CPFEM_ffn(:,:,g,i,e)
 subFrac = 0.0_pReal
 subStep = 1.0_pReal
 nCutbacks = 0_pInt
 maxCutbacks = 0_pInt
 Fg_current = CPFEM_ffn(:,:,g,i,e)                                 ! initialize to start of inc
 Fp_current = CPFEM_Fp_old(:,:,g,i,e)
 call math_invert3x3(Fp_current,inv,det,error)
 Fe_current = math_mul33x33(Fg_current,inv)
! state_current = state_new
 success = .false.                                                 ! pretend cutback
 dt_aim = 0.0_pReal                                                ! prevent initial Lp interpolation
 Temperature=CPFEM_Temperature(i,e)
! 
! begin the cutback loop
 do while (subStep > subStepMin)                                   ! continue until finished or too much cut backing
   if (success) then                                               ! wind forward
     Fg_current = Fg_aim
     Fe_current = Fe_new
     Fp_current = Fp_new
!     state_current = state_new
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
!     write (6,'(a,/,3(4(f9.3,x)/))') 'state current / MPa',state_current/1e6_pReal
     write (6,'(a,/,3(4(f9.3,x)/))') 'state new / MPa',state_new/1e6_pReal
     write (6,'(a,/,3(3(f12.7,x)/))') 'Fe current',Fe_current(1:3,:)
     write (6,'(a,/,3(3(f12.7,x)/))') 'Fp current',Fp_current(1:3,:)
     write (6,'(a,/,3(3(f12.7,x)/))') 'Lp (old=new guess)',Lp(1:3,:)
     write (6,'(a20,f,x,a2,x,f)') 'integrating from ',subFrac,'to',(subFrac+subStep)
!$OMP END CRITICAL (write2out)
   endif
!
   call TimeIntegration(msg,Lp,Fp_new,Fe_new,Tstar_v,P,state_new,dt_aim,e,i,g,Temperature,Fg_aim,Fp_current)
!

   if (msg == 'ok') then
     subFrac = subFrac + subStep
     subStep = min(1.0_pReal-subFrac, subStep*2.0_pReal) ! accelerate
     nCutbacks = 0_pInt                   ! reset cutback counter
     success = .true.                     ! keep current Lp
   else
     nCutbacks = nCutbacks + 1            ! record additional cutback
     maxCutbacks = max(nCutbacks,maxCutbacks)! remember maximum number of cutbacks
     subStep = subStep / 2.0_pReal        ! cut time step in half
     success = .false.                    ! force Lp interpolation
!     if (debugger) then
!$OMP CRITICAL (write2out)
       write (6,*) '>>>>>>>>>>>>>>>>>>>> cutback <<<<<<<<<<<<<<<<<<<<<<'
       write (6,*) 'Element, Ip:', e, i
       write (6,*) msg
!$OMP END CRITICAL (write2out)
!     endif
!
   endif
 enddo  ! potential substepping
!
!$OMP CRITICAL (cutback)
 debug_cutbackDistribution(min(nCutback,maxCutbacks)+1) = debug_cutbackDistribution(min(nCutback,maxCutbacks)+1)+1
!$OMP END CRITICAL (cutback)
!
! debugger = .false.
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
     P,&                ! 1nd PK stress (taken as initial guess if /= 0)
     state,&            ! current microstructure at end of time inc (taken as guess if /= 0)
     dt,&               ! time increment
     cp_en,&            ! element number
     ip,&               ! integration point number
     grain,&            ! grain number
     Temperature,&      ! temperature
     Fg_new,&           ! new total def gradient
     Fp_old)            ! former plastic def gradient
!     state_current)         ! former microstructure
 use prec
 use debug
 use mesh, only: mesh_element
 use constitutive, only: constitutive_Nstatevars,constitutive_Microstructure,&
                         constitutive_homogenizedC,constitutive_LpAndItsTangent
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
 real(pReal), dimension(3,3) :: P !,Tstar
 real(pReal), dimension(3,3) :: Lp,Lpguess,Lpguess_old,Rinner,Rinner_old,A,B,BT,AB,BTA
 real(pReal), dimension(3,3,3,3) :: C
 real(pReal), dimension(constitutive_Nstatevars(grain, ip, cp_en)) :: state
!
 msg = 'ok'  ! error-free so far
 eye2 = math_identity2nd(9)

 call math_invert3x3(Fp_old,invFp_old,det,failed) ! inversion of Fp_old
 if (failed) then
    msg = 'inversion Fp_old'
    return
 endif

 A = math_mul33x33(transpose(invFp_old), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_old)))
!
! if (all(state == 0.0_pReal)) state = state_current ! former state guessed, if none specified
! iOuter = 0_pInt                                    ! outer counter
!
 if (debugger) then
!$OMP CRITICAL (write2out)
   write (6,'(a,/,3(3(f12.7,x)/))') 'Fg to be calculated',Fg_new
!$OMP END CRITICAL (write2out)
 endif
!
!Outer: do                ! outer iteration: State
!         iOuter = iOuter+1
!         if (debugger) then
!!$OMP CRITICAL (write2out)
!           write (6,'(a,i3)') '---outer ',iOuter
!           write (6,'(a,/,3(4(f9.3,x)/))') 'state old / MPa',state_old/1e6_pReal
!           write (6,'(a,/,3(4(f9.3,x)/))') 'state / MPa',state/1e6_pReal
!           write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
!!$OMP END CRITICAL (write2out)
!         endif
!
!         if (iOuter > nOuter) then
!           msg = 'limit Outer iteration'
!!$OMP CRITICAL (out)
!           debug_OuterLoopDistribution(nOuter) = debug_OuterLoopDistribution(nOuter)+1
!!$OMP END CRITICAL (out)
!           return
!         endif
         call constitutive_Microstructure(state,Temperature,grain,ip,cp_en)
         C_66 = constitutive_HomogenizedC(state, grain, ip, cp_en)
         C = math_Mandel66to3333(C_66)       ! 4th rank elasticity tensor
!
         iInner = 0_pInt
         leapfrog = 1.0_pReal                ! correction as suggested by invdRdLp-step
         maxleap = 1024.0_pReal              ! preassign maximum acceleration level
!
		 Lpguess_old = Lpguess               ! consider present Lpguess good
!
Inner: do              ! inner iteration: Lp
         iInner = iInner+1
!         if (debugger) then
!!$OMP CRITICAL (write2out)
!           write (6,'(a,i3)') 'inner ',iInner
!		   if (iInner < 3) then
!             write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
!		   endif
!!$OMP END CRITICAL (write2out)
!         endif
         if (iInner > nInner) then           ! too many loops required
           Lpguess = Lpguess_old             ! do not trust the last update but resort to former one
           msg = 'limit Inner iteration'
!$OMP CRITICAL (in)
           debug_InnerLoopDistribution(nInner) = debug_InnerLoopDistribution(nInner)+1
!$OMP END CRITICAL (in)
           return
         endif
!
         B = math_i3 - dt*Lpguess
         BT = transpose(B)
         AB = math_mul33x33(A,B)
         BTA = math_mul33x33(BT,A)
         Tstar_v = 0.5_pReal*math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB)-math_I3))
!         Tstar = math_Mandel6to33(Tstar_v)
         p_hydro=(Tstar_v(1)+Tstar_v(2)+Tstar_v(3))/3.0_pReal
         forall(i=1:3) Tstar_v(i) = Tstar_v(i)-p_hydro                ! subtract hydrostatic pressure
         call constitutive_LpAndItsTangent(Lp,dLp, &
                                           Tstar_v,state,Temperature,grain,ip,cp_en)
!			    
         Rinner = Lpguess - Lp                                        ! update current residuum
!
         if (.not.(any(Rinner/=Rinner)) .and. &                       ! exclude any NaN in residuum
             ( (maxval(abs(Rinner)) < abstol_Inner) .or. &            ! below abs tol .or.
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
                 write (6,'(a,/,3(4(f9.3,x)/))') 'state / MPa',state/1e6_pReal
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
           if (iInner > 1 .and. leapfrog < maxleap) &
		     leapfrog = 2.0_pReal * leapfrog                          ! accelerate if ok
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
!       ROuter = state - state_old - &
!                dt*constitutive_dotState(Tstar_v,state,Temperature,&
!                                         grain,ip,cp_en)          ! residuum from evolution of microstructure
!       state = state - ROuter                                           ! update of microstructure
!
!	   if (iOuter==nOuter) then
!!$OMP CRITICAL (write2out)
!	      write (6,*) 'Terminated outer loop at el,ip,grain',cp_en,ip,grain
!!$OMP END CRITICAL (write2out)
!	      exit Outer
!	   endif
!       if (maxval(abs(Router/state),state /= 0.0_pReal) < reltol_Outer) exit Outer
!     enddo Outer
!
!!$OMP CRITICAL (out)
! debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!!$OMP END CRITICAL (out)
 invFp_new = math_mul33x33(invFp_old,B)
 call math_invert3x3(invFp_new,Fp_new,det,failed)
 if (failed) then
   msg = 'inversion Fp_new^-1'
   return
 endif
!
! if (wantsConstitutiveResults) then     ! get the post_results upon request
!   results = 0.0_pReal
!   results = constitutive_post_results(Tstar_v,state,Temperature,dt,grain,ip,cp_en)
! endif
!
 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal)     ! regularize Fp by det = det(InvFp_new) !!
 forall (i=1:3) Tstar_v(i) = Tstar_v(i)+p_hydro ! add hydrostatic component back
 Fe_new = math_mul33x33(Fg_new,invFp_new)              ! calc resulting Fe
! P = math_mul33x33(Fe_new,math_mul33x33(Tstar,transpose(invFp_new)))    ! first PK stress
 P = math_mul33x33(Fe_new,math_mul33x33(math_Mandel6to33(Tstar_v),transpose(invFp_new)))    ! first PK stress

 return
!
 END SUBROUTINE
!
 END MODULE
!##############################################################

