!    ---------------------------
 MODULE CPFEM     
!    ---------------------------
!    *** CPFEM engine ***
!
 use prec, only: pReal,pInt
 implicit none
!
!    ****************************************************************
!    *** General variables for the material behaviour calculation ***
!    ****************************************************************  
 real(pReal), allocatable :: CPFEM_stress_all       (:,:,:)
 real(pReal), allocatable :: CPFEM_jacobi_all       (:,:,:,:)
 real(pReal), allocatable :: CPFEM_ffn_all          (:,:,:,:)
 real(pReal), allocatable :: CPFEM_ffn1_all         (:,:,:,:)
 real(pReal), allocatable :: CPFEM_results          (:,:,:,:)
 real(pReal), allocatable :: CPFEM_ini_ori          (:,:,:,:)
 real(pReal), allocatable :: CPFEM_sigma_old        (:,:,:,:)
 real(pReal), allocatable :: CPFEM_sigma_new        (:,:,:,:)
 real(pReal), allocatable :: CPFEM_Fp_old           (:,:,:,:,:)
 real(pReal), allocatable :: CPFEM_Fp_new           (:,:,:,:,:)
 real(pReal), allocatable :: constitutive_state_old (:,:,:,:)
 real(pReal), allocatable :: constitutive_state_new (:,:,:,:)
 real(pReal), allocatable :: CPFEM_g_old            (:,:,:,:)
 real(pReal), allocatable :: CPFEM_g_new            (:,:,:,:)
 real(pReal), allocatable :: CPFEM_jaco_old         (:,:,:,:)
 real(pReal), allocatable :: CPFEM_mat              (:,:)
 integer(pInt) :: CPFEM_inc_old    = 0_pInt
 integer(pInt) :: CPFEM_subinc_old = 1_pInt
 integer(pInt) :: CPFEM_first_call = 1_pInt
 integer(pInt) :: CPFEM_Nresults   = 4_pInt

 CONTAINS

!***********************************************************************
!***    This routine checks for initialization, variables update and ***
!***    calls the actual material model                              ***
!***********************************************************************
 subroutine cpfem_general(ffn, ffn1, ndi, CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_dt, CPFEM_en, CPFEM_in)
!
 use prec, only: pReal,pInt
 use CPFEM, only: CPFEM_ffn_all, CPFEM_ffn1_all, CPFEM_inc_old
 use IO, only: IO_error
 implicit none
!
 real(pReal)   ffn(3,3), ffn1(3,3), CPFEM_dt
 integer(pInt) ndi, CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_en, CPFEM_in
!
! initialization step
 if (CPFEM_first_call==1_pInt) then
! three dimensional stress state ?
    if (CPFEM_ndi/=3_pInt) then
        call IO_error(300)
    endif
    call IO_allocation()
    call mesh_allocation()
    call constitutive_allocation()
    call math_allocation()
    call IO_allocation()
    call CPFEM_allocation()
    CPFEM_first_call=0_pInt
 endif
! not a new increment
 if (CPFEM_inc==CPFEM_inc_old) then
! case of a new subincrement:update starting with subinc 2
     if (CPFEM_subinc > CPFEM_subinc_old) then
        CPFEM_sigma_old        = CPFEM_sigma_new
        CPFEM_Fp_old           = CPFEM_Fp_new
        constitutive_state_old = constitutive_state_new
        CPFEM_g_old            = CPFEM_g_new
        CPFEM_subinc_old       = CPFEM_subinc
    endif
    return
! case of a new increment
 else
    CPFEM_sigma_old         = CPFEM_sigma_new
    CPFEM_Fp_old            = CPFEM_Fp_new
    constitutive_state_old  = constitutive_state_new
    CPFEM_g_old             = CPFEM_g_new
    CPFEM_inc_old           = CPFEM_inc
    CPFEM_subinc_old        = 1_pInt
    CPFEM_timefactor_max    = 0.0_pReal
 endif
!
! get cp element number for fe element number
 cp_en=mesh_??(CPFEM_en)!ÄÄÄ
 CPFEM_ffn_all(:,:,CPFEM_in, cp_en)  = ffn
 CPFEM_ffn1_all(:,:,CPFEM_in, cp_en) = ffn1
 call CPFEM_general_material(CPFEM_cn, CPFEM_dt, cp_en, CPFEM_in)
 return
 end


!***********************************************************************
!***    This routine allocates the arrays defined in module CPFEM    ***
!***    and initializes them                                         ***
!***********************************************************************
 subroutine CPFEM_allocation()
!
 use prec, only: pReal,pInt
 use IO, only: IO_error
 use math
 use mesh
 use constitutive
!
 implicit none
!
 integer(pInt) i
!
!    *** mpie.marc parameters ***
 allocate(CPFEM_ffn_all(3,3,mesh_Nips,mesh_Nelems))
 allocate(CPFEM_ffn1_all(3,3,mesh_Nips,mesh_Nelems))
 allocate(CPFEM_stress_all(6,mesh_Nips,mesh_Nelems))
 allocate(CPFEM_jacobi_all(6,6,mesh_Nips,mesh_Nelems))
 CPFEM_ffn_all    = 0.0_pReal
 CPFEM_ffn1_all   = 0.0_pReal
 forall(i=1:3)
    CPFEM_ffn_all(i,i,:,:,:)  = 1.0_pReal 
    CPFEM_ffn1_all(i,i,:,:,:) = 1.0_pReal
 endforall 
 CPFEM_stress_all = 0.0_pReal
 CPFEM_jacobi_all = 0.0_pReal
!
!    *** User defined results ***
 allocate(CPFEM_results(CPFEM_Nresults,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 CPFEM_results = 0.0_pReal
!
!    *** Initial orientations ***
! allocate(CPFEM_ini_ori(3,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
! CPFEM_ini_ori = 0.0_pReal
!
!    *** Second Piola-Kirchoff stress tensor at (t=t0) and (t=t1) ***
 allocate(CPFEM_sigma_old(6,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 allocate(CPFEM_sigma_new(6,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 CPFEM_sigma_old = 0.0_pReal
 CPFEM_sigma_new = 0.0_pReal
!
!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***  
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 CPFEM_Fp_old = 0.0_pReal 
 CPFEM_Fp_new = 0.0_pReal 
 forall(i=1:3)
    CPFEM_Fp_old(i,i,:,:,:) = 1.0_pReal 
    CPFEM_Fp_new(i,i,:,:,:) = 1.0_pReal
 endforall 
!    
! QUESTION: would it be wise to outsource these to _constitutive_ ??
!    *** Slip resistances at (t=t0) and (t=t1) ***  
 allocate(constitutive_state_old(constitutive_Nstatevars,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 allocate(constitutive_state_new(constitutive_Nstatevars,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 state_tauc_slip_old = 0.0_pReal
 state_tauc_slip_new = 0.0_pReal

!    *** Cumulative shear at (t=t0) and (t=t1) *** 
! QUESTION which nslip to use here ?!?
 allocate(CPFEM_g_old(constitutive_maxNslip,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))  
 allocate(CPFEM_g_new(constitutive_maxNslip,constitutive_maxNgrains,mesh_Nips,mesh_Nelems))
 CPFEM_g_old = 0.0_pReal
 CPFEM_g_new = 0.0_pReal
!
!    *** Old jacobian (consistent tangent) ***
 allocate(CPFEM_jaco_old(6,6,mesh_Nips,mesh_Nelems))
 CPFEM_jaco_old = 0.0_pReal
!
!    *** Output to MARC output file ***
 write(6,*)
 write(6,*) 'Arrays allocated:'
 write(6,*) 'CPFEM_ffn_all:       ', shape(CPFEM_ffn_all)
 write(6,*) 'CPFEM_ffn1_all:      ', shape(CPFEM_ffn1_all)
 write(6,*) 'CPFEM_stress_all:    ', shape(CPFEM_stress_all)
 write(6,*) 'CPFEM_jacobi_all:    ', shape(CPFEM_jacobi_all)
 write(6,*) 'CPFEM_results:       ', shape(CPFEM_results)
 write(6,*) 'CPFEM_thickness:     ', shape(CPFEM_thickness)
! write(6,*) 'CPFEM_ini_ori:       ', shape(CPFEM_ini_ori)
 write(6,*) 'CPFEM_sigma_old:     ', shape(CPFEM_sigma_old)
 write(6,*) 'CPFEM_sigma_new:     ', shape(CPFEM_sigma_new)
 write(6,*) 'CPFEM_Fp_old:        ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:        ', shape(CPFEM_Fp_new)
 write(6,*) 'constitutive_state_old: ', shape(constitutive_state_old)
 write(6,*) 'constitutive_state_new: ', shape(constitutive_state_new)
 write(6,*) 'CPFEM_g_old:         ', shape(CPFEM_g_old)
 write(6,*) 'CPFEM_g_new:         ', shape(CPFEM_g_new)
 write(6,*) 'CPFEM_jaco_old:      ', shape(CPFEM_jaco_old)
 write(6,*)
 call flush(6)
 return
 end     
!
!
 subroutine CPFEM_general_material(&
     CPFEM_cn,&       ! Cycle number
     CPFEM_dt,&       ! Time increment (dt)
     CPFEM_en,&       ! Element number
     CPFEM_in)        ! Integration point number
!***********************************************************************
!***    This routine calculates the material behaviour    ***
!***********************************************************************
 use prec, only: pReal,pInt
 use IO, only: IO_error
 use math
 use mesh
 use constitutive
!
 implicit none
!
!    *** Definition of variables ***
 integer(pInt) CPFEM_cn, CPFEM_en ,CPFEM_in
 real(pReal) CPFEM_dt, CPFEM_s(6), CPFEM_d(6, 6), CPFEM_ffn(3,3),CPFEM_ffn1(3,3)
! QUESTION which nslip to use?
 real(pReal) Fp_old(3,3), tauc_slip_old(constitutive_maxNslip), tauc_slip_new(constitutive_maxNslip), g_old(constitutive_maxNslip)
 real(pReal) g_new(constitutive_maxNslip), Tstar_v(6), Fp_new(3,3), cs(6), phi1mis(2), PHImis(2), phi2mis(2), cd(6,6)
 real(pReal) ori_mat(3,3),hh6(6,6)
 integer(pInt) jpara,nori
 real(pReal) phi1, PHI, phi2, scatter, vf, alpha1, alpha2, beta1, &
             beta2, phi1_s, PHI_s, phi2_s, p10, P0, p20, p11, P1, p21, &
             dgmax,dgmaxc , orimis
 integer(pInt) i, iori, iconv, ising, icut
!    *** Numerical parameters ***
!    *** How often the jacobian is recalculated ***
 integer (pInt), parameter :: ijaco = 5_pInt
!    *** Reference shear rate for the calculation of CPFEM_timefactor ***
 real (pReal), parameter :: dgs = 0.01_pReal
!
!    *** Flag for recalculation of jacobian ***
 jpara = 1_pInt      
! get cp element number for fe element number
 cp_en=mesh_??(CPFEM_en)!ÄÄÄ
! get number of grains from cp element number and integration point number
 nori = constitutive_???(cp_en, CPFEM_in)  !ÄÄÄ
! 
!
 CPFEM_s=0
 CPFEM_d=0
!
!    *** Loop over all the components ***
 do iori=1,nori
!
!    *** Initialization of the matrices for t=t0 ***
!    Fp_old        = CPFEM_Fp_old(:,:,iori,CPFEM_in,cp_en)
!    tauc_slip_old = constitutive_state_old(:,iori,CPFEM_in,cp_en)
!    tauc_slip_new = tauc_slip_old
!    g_old         = CPFEM_g_old(:,iori,CPFEM_in,cp_en)
!    Tstar_v       = CPFEM_sigma_old(:,iori,CPFEM_in,cp_en)
! data from constitutive?
    vf            = constitutive_vol(iori,CPFEM_in,cp_en) !ÄÄÄ

!    *** Calculation of the solution at t=t1 ***
    if (modulo(CPFEM_cn,ijaco)==0) then !ÄÄÄ
        call CPFEM_stress(cs, cd, CPFEM_dt,cp_en,CPFEM_in, iori, ising, icut, iconv, dgmaxc, 1_pInt)
!
!
!        call CPFEM_stress(CPFEM_dt,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new,
!     &     g_old,g_new,tauc_slip_old, 
!     &     tauc_slip_new,
!     &     Tstar_v,cs,cd,p11,P1,p21,dgmaxc,1,iconv,ising,
!     &     icut,CPFEM_en,CPFEM_in,CPFEM_inc)
!    ***   Evaluation of ising      ***
!    *** ising=2 => singular matrix in jacobi calculation ***
!    ***      => use old jacobi        ***
        if (ising==2) jpara=0
!    *** Calculation of the consistent tangent ***
        CPFEM_d=CPFEM_d+vf*cd
    else
       call CPFEM_stress(cs, cd, CPFEM_dt,cp_en,CPFEM_in, iori, ising, icut, iconv, dgmaxc, 0_pInt)
!    call CPFEM_stress(CPFEM_tinc,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new, 
!     &     g_old,g_new,tauc_slip_old,
!     &     tauc_slip_new,
!     &     Tstar_v,cs,hh6,p11,P1,p21,dgmaxc,0,iconv,
!     &     ising,icut,CPFEM_en,CPFEM_in,CPFEM_inc)
       jpara=0
    endif
!    *** Cases of unsuccessful calculations *** 
!    ***   Evaluation od ising  ***
!    *** ising!=0 => singular matrix ***
    if (ising==1) then
        write(6,*) 'Singular matrix!'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(700)
        CPFEM_timefactor=1.e5_pReal
        return
    endif
!    ***   Evaluation of icut   ***
!    *** icut!=0 => too many cutbacks ***
    if (icut==1) then
        write(6,*) 'Too many cutbacks'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(600)
        CPFEM_timefactor=1.e5_pReal
        return
    endif
!    ***  Evaluation of iconv ***
!    *** iconv!=0 => no convergence ***
    if (iconv==1) then
        write(6,*) 'Inner loop did not converged!'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(600)
        CPFEM_timefactor=1.e5_pReal
        return
    else if (iconv==2) then
        write(6,*) 'Outer loop did not converged!'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(600)
        CPFEM_timefactor=1.e5_pReal
        return
    endif
!    *** Update the differents matrices for t=t1 ***
!    CPFEM_Fp_new(:,:,iori,CPFEM_in,cp_en)      = Fp_new
!    state_tauc_slip_new(:,iori,CPFEM_in,cp_en) = tauc_slip_new
!    CPFEM_g_new(:,iori,CPFEM_in,cp_en)         = g_new
!    CPFEM_sigma_new(:,iori,CPFEM_in,cp_en)     = Tstar_v
!
!    *** Calculation of the misorientation ***    
!phi1mis(1)=p10
! PHImis(1)=P0
! phi2mis(1)=p20
! phi1mis(2)=p11
! PHImis(2)=P1
! phi2mis(2)=p21
! call CPFEM_misori(phi1mis,PHImis,phi2mis,orimis)
!
!    *** Update the results plotted in MENTAT ***
! CPFEM_results(1,iori,cp_en,CPFEM_in) = p11
! CPFEM_results(2,iori,cp_en,CPFEM_in) = P1
! CPFEM_results(3,iori,cp_en,CPFEM_in) = p21
! CPFEM_results(4,iori,cp_en,CPFEM_in) = sum(g_new)
!
!    *** Evaluation of the maximum shear ***
    dgmax=max(dgmax,dgmaxc)
!    *** Evaluation of the average Cauchy stress ***    
    CPFEM_s=CPFEM_s+vf*cs
 enddo
!    *** End of the loop over the components ***
!    *************************************
!    *** End of the CP-FEM Calculation *** 
!    *************************************
!    *** Restoration of the old jacobian if necessary ***
 if (jpara==0) then
    CPFEM_d=CPFEM_jaco_old(:,:,CPFEM_in,cp_en)
 else
!    *** Store the new jacobian ***     
    CPFEM_jaco_old(:,:,CPFEM_in,cp_en)=CPFEM_d
 endif    
!    *** Calculate timefactor ***
 CPFEM_timefactor=dgmax/dgs
!
 return
 end


!call CPFEM_stress(cs, cd, CPFEM_dt,cp_en,CPFEM_in, ising, icut, iconv, dgmaxc, 1)
 subroutine CPFEM_stress(&
     cs,&               ! stress vector
     cd,&               ! Jacoby matrix
     CPFEM_dt,&         ! Time increment (dt)
     cp_en,&            ! Element number
     CPFEM_in,&         ! Integration point number
     iori,&             ! number of orintation
     ising,&            ! flag for singular matrix
     icut,&             ! flag for too many cut backs
     iconv,&            ! flag for non convergence
     dgmaxc,&           ! maximum shear
     isjaco)            ! flag whether to calculate Jacoby matrix
!********************************************************************
! This routine calculates the stress for a single component
! and manages the independent time incrmentation
!********************************************************************
 use prec, only: pReal,pInt
 use CPFEM, only: CPFEM_ffn_all, CPFEM_ffn1_all
 implicit none
!
!    *** Definition of variables ***
 integer(pInt) isjaco,iconv,ising,icut,CPFEM_en,CPFEM_in,CPFEM_inc
 real(pReal)   CPFEM_tinc,CPFEM_ffn(3,3),CPFEM_ffn1(3,3),Fp_old(3,3)
 real(pReal)   Fp_new(3,3),g_old(nslip),g_new(nslip)
 real(pReal)   tauc_slip_old(nslip),tauc_slip_new(nslip)
 real(pReal)   Tstar_v(6)
 real(pReal)   cs(6),dcs_de(6,6),phi1,PHI,phi2,dgmaxc
 integer(pInt) jcut
 real(pReal)   Tstar_v_h(6),tauc_slip_new_h(nslip)
 real(pReal)   dt_i,delta_Fg(3,3),Fg_i(3,3)
 real(pReal)   tauc_slip_new_i(nslip),time,mm(6,6)

!    *** Numerical parameters ***
 integer(pInt), parameter :: ncut=7_pInt
 icut=0
!
!    *** Initialization of the matrices for t=t0 ***
    Fp_old     = CPFEM_Fp_old(:,:,iori,CPFEM_in,cp_en)
    Fp_new     = 0_pReal
    state_old  = constitutive_state_old(:,iori,CPFEM_in,cp_en)
    state_new  = state_old
    g_old      = CPFEM_g_old(:,iori,CPFEM_in,cp_en)
    g_new      = 0_pReal
    Tstar_v    = CPFEM_sigma_old(:,iori,CPFEM_in,cp_en)
    CPFEM_ffn  = CPFEM_ffn_all(:,:,CPFEM_in,cp_en)
    CPFEM_ffn1 = CPFEM_ffn1_all(:,:,CPFEM_in,cp_en)
!
!    *** First attempt to calculate Tstar and tauc with initial timestep ***
 Tstar_v_h=Tstar_v
 state_new_h=state_new
 call CPFEM_stress_int(cs, cd, CPFEM_dt, cp_en,CPFEM_in, ising, icut, iconv, dgmaxc, isjaco, phi1, PHI, phi2,&
                       CPFEM_ffn, CPFEM_ffn1,Fp_old,Fp_new,g_old,g_new,state_old, state_new, Tstar_v)
 if ((iconv==0).AND.(ising==0)) then
!    *** Update the differents matrices for t=t1 ***
    CPFEM_Fp_new(:,:,iori,CPFEM_in,cp_en)          = Fp_new
    constituitive_state_new(:,iori,CPFEM_in,cp_en) = state_new
    CPFEM_g_new(:,iori,CPFEM_in,cp_en)             = g_new
    CPFEM_sigma_new(:,iori,CPFEM_in,cp_en)         = Tstar_v
!    *** Update the results plotted in MENTAT ***
    CPFEM_results(1,iori,CPFEM_in,cp_en) = phi1
    CPFEM_results(2,iori,CPFEM_in,cp_en) = PHI
    CPFEM_results(3,iori,CPFEM_in,cp_en) = phi2
    CPFEM_results(4,iori,CPFEM_in,cp_en) = sum(g_new)
    return
 endif
!
!    *** Calculation of stress and resistences with a cut timestep ***
!    ***      when first try did not converge    ***  
 jcut=1_pInt
 dt_i=0.5_pReal*CPFEM_dt
 delta_Fg=0.5_pReal*(CPFEM_ffn1-CPFEM_ffn)
 Fg_i=CPFEM_ffn+delta_Fg
 Tstar_v=Tstar_v_h
 state_new_i=state_new_h
!    *** Start time ***
 time=dt_i
 do while (time<=CPFEM_dt)
    call CPFEM_stress_int(cs, cd, dt_i, cp_en,CPFEM_in, ising, icut, iconv, dgmaxc, isjaco, phi1, PHI, phi2,&
                       CPFEM_ffn, Fg_i,Fp_old,Fp_new,g_old,g_new,state_old, state_new_i, Tstar_v)
!    call CPFEM_stress_int(time,CPFEM_ffn,Fg_i,Fp_old,Fp_new,g_old,
!     &    g_new,tauc_slip_old,
!     &    tauc_slip_new_i,
!     &    Tstar_v,cs,mm,phi1,PHI,
!     &    phi2,dgmaxc,0_pInt,iconv,ising,CPFEM_en,
!     &    CPFEM_in,CPFEM_inc)
    if ((iconv==0).AND.(ising==0)) then
        time=time+dt_i
        Fg_i=Fg_i+delta_Fg
        Tstar_v_h=Tstar_v
        state_new_h=state_new_i
    else
        jcut=jcut+1_pInt
        if (jcut.GT.ncut) then
            icut=1_pInt
            return
        endif
        dt_i=0.5_pReal*dt_i
        time=time-dt_i
        delta_Fg=0.5_pReal*delta_Fg
        Fg_i=Fg_i-delta_Fg
        Tstar_v=Tstar_v_h
        tauc_slip_new_i=tauc_slip_new_h
    endif
 enddo
!
!    *** Final calculation of stress and resistences withb full timestep ***
 state_new=state_new_i
 call CPFEM_stress_int(cs, cd, CPFEM_dt, cp_en,CPFEM_in, ising, icut, iconv, dgmaxc, isjaco, phi1, PHI, phi2,&
                       CPFEM_ffn, CPFEM_ffn1,Fp_old,Fp_new,g_old,g_new,state_old, state_new, Tstar_v)
! call CPFEM_stress_int(CPFEM_tinc,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new,
!     &      g_old,g_new,tauc_slip_old,
!     &      tauc_slip_new,
!     &      Tstar_v,cs,dcs_de,phi1,PHI,phi2,dgmaxc,
!     &      isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc) 
!    *** Update the differents matrices for t=t1 ***
 CPFEM_Fp_new(:,:,iori,CPFEM_in,cp_en)          = Fp_new
 constituitive_state_new(:,iori,CPFEM_in,cp_en) = state_new
 CPFEM_g_new(:,iori,CPFEM_in,cp_en)             = g_new
 CPFEM_sigma_new(:,iori,CPFEM_in,cp_en)         = Tstar_v
!    *** Update the results plotted in MENTAT ***
 CPFEM_results(1,iori,CPFEM_in,cp_en) = phi1
 CPFEM_results(2,iori,CPFEM_in,cp_en) = PHI
 CPFEM_results(3,iori,CPFEM_in,cp_en) = phi2
 CPFEM_results(4,iori,CPFEM_in,cp_en) = sum(g_new)
 return
 end


! call CPFEM_stress_int(cs, cd, CPFEM_dt, cp_en,CPFEM_in, ising, icut, iconv, dgmaxc, isjaco,&
!                       CPFEM_ffn, CPFEM_ffn1,Fp_old,Fp_new,g_old,g_new,state_old, state_new, Tstar_v)

 subroutine CPFEM_stress_int(&
     cs,&         ! Cauchy stress vector
     dcs_de,&     ! Consistent tangent
     dt,&         ! Time increment
     cp_en,&      ! cp element number
     CPFEM_in,&   ! integration point number
     ising,&      ! flag for singular matrix
     icut,&       ! flag for too many cut backs
     iconv,&      ! flag for non convergence
     dgmaxc,&     ! maximum shear
     isjaco,&     ! flag whether to calculate Jacoby matrix
     phi1,&       ! Euler angle
     PHI,&        ! Euler angle
     phi2,&       ! Euler angle
     Fg_old,&     ! Old global deformation gradient
     Fg_new,&     ! New global deformation gradient
     Fp_old,&     ! Old plastic deformation gradient
     Fp_new,&     ! New plastic deformation gradient
     g_old,&      ! Old cumulative plastic strain of a slip system
     g_new,&      ! New cumulative plastic strain of a slip system
     state_old,&  ! Old resistence of a slip system
     state_new,&  ! New resistence of a slip system
     Tstar_v)     ! Second Piola-Kirschoff stress tensor
!********************************************************************
! This routine calculates the stress for a single component
! it is based on the paper by Kalidindi et al.:
! J. Mech. Phys, Solids Vol. 40, No. 3, pp. 537-569, 1992
! it is modified to use anisotropic elasticity matrix
!********************************************************************
 use prec
 implicit none

!    *** Definition of variables *** 
 integer(pInt) isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc
 real(pReal) dt,Fg_old(3,3),Fg_new(3,3),Fp_old(3,3),Fp_new(3,3),
     &      g_old(nslip),g_new(nslip),
     &      tauc_slip_old(nslip),tauc_slip_new(nslip),
     &      Tstar_v(6),
     &      cs(6),dcs_de(6,6),phi1,PHI,phi2,dgmaxc
 integer(pInt) ic     
 real(pReal) gdot_slip(nslip),Fe(3,3),R(3,3),
     &      U(3,3),de(3,3),tauc2(nslip),Fp2(3,3),
     &      sgm2(6),cs1(6),dF(3,3),Fg2(3,3),dev(6)    
!    *** Numerical parameters ***
 real(pReal), parameter :: pert_ct=1.0e-5_pReal  
! maximum shear rate 
 dgmaxc = 0
!    *** Error treatment ***
 iconv  = 0
 ising  = 0   

!    *********************************************
!    ***  Calculation of the new Cauchy stress ***
!    *********************************************

!    *** Call Newton-Raphson method ***
 call NEWTON_RAPHSON(dt,Fg_old,Fg_new,Fp_old,Fp_new,Fe,gdot,state_old,state_new,Tstar_v,cs,iconv,ising)
! 
!    *** Calculation of the new orientation ***
 call math_pDecomposition(Fe,U,R,ising)
 if (ising==1) then
    return
 endif
 call math_RtoEuler(transpose(R),phi1,PHI,phi2)
!
!    *** Evaluation of the maximum slip shear ***
 dgmaxc=maxval(abs(gdot_slip*dt))
 g_new=g_old+abs(gdot_slip)*dt
!
!    *** Choice of the calculation of the consistent tangent ***     
 if (isjaco==0) return
!
!    *********************************************
!    *** Calculation of the consistent tangent ***
!    *********************************************
!
!    *** Calculation of the consistent tangent with perturbation ***
!    *** Perturbation on the component of Fg ***    
 do ic=1,6    
!
!    *** Method of small perturbation
    dev=0
    if(ic<=3) dev(ic) = pert_ct
    if(ic>3)  dev(ic) = pert_ct/2
    call math_conv6to33(dev,de)
    dF=matmul(de,Fg_old) 
    Fg2=Fg_new+dF
    sgm2=Tstar_v
    state2=state_new

!    *** Calculation of the perturbated Cauchy stress ***
    call NEWTON_RAPHSON(dt,Fg_old,Fg2,Fp_old,Fp2,Fe,gdot,state_old,tauc2,sgm2,cs1,iconv,ising)
!  
!    *** Consistent tangent ***
    dcs_de(:,ic)=(cs1-cs)/pert_ct     
 enddo
!
 return
 end
!
!
 subroutine NEWTON_RAPHSON(
     &dt, 
     &Fg_old,     
     &Fg_new,     
     &Fp_old,     
     &Fp_new,
     &Fe,       
     &gdot_slip,    
     &state_old,  
     &state_new,  
     &Tstar_v, 
     &cs,
     &iconv,
     &ising          
     &) 
!***********************************************************************
!***        NEWTON-RAPHSON Calculation                               ***
!***********************************************************************
 use prec
 implicit none

!    *** Definition of variables *** 
 integer(pInt) isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc
 real(pReal) dt,Fg_old(3,3),Fg_new(3,3),Fp_old(3,3),Fp_new(3,3),
     &      g_old(nslip),g_new(nslip),
     &      tauc_slip_old(nslip),tauc_slip_new(nslip),
     &      Tstar_v(6),cs(6),dcs_de(6,6),phi1,PHI,phi2,dgmaxc
 integer(pInt) i,j,k,iouter,iinner,ijac,ic      
 real(pReal) invFp_old(3,3),det,A(3,3),Estar0_v(6),Tstar0_v(6),
     &      mm(3,3),mm1(3,3),vv(6),Dslip(6,nslip),
     &      tau_slip(nslip),gdot_slip(nslip),
     &      R1(6),norm1,Tstar_v_per(6),R1_per(6),
     &      Jacobi(6,6),invJacobi(6,6),dTstar_v(6),R2(nslip),
     &      dtauc_slip(nslip),norm2,dLp(3,3),
     &      Estar(3,3),Estar_v(6),invFp_new(3,3),
     &      invFp2(3,3),Lp(3,3),Fe(3,3),
     &      R(3,3),U(3,3),dgdot_dtaucslip(nslip)
 real(pReal) de(3,3),dev(6),tauc2(nslip),fp2(3,3),
     &      sgm2(6),cs1(6),df(3,3),
     &      fg2(3,3),tauc_old(nslip),crite,tol_in,tol_out

!    *** Numerical parameters ***
 integer(pInt), parameter :: nouter    = 50_pInt
 real(pReal),   parameter :: tol_outer = 1.0e-4_pReal
 integer(pInt), parameter :: ninner    = 2000_pInt
 real(pReal),   parameter :: tol_inner = 1.0e-3_pReal   
 real(pReal),   parameter :: eta       = 13.7_pReal 
 integer(pInt), parameter :: numerical = 0_pInt
 real(pReal),   parameter :: pert_nr   = 1.0e-8_pReal  
 crite=eta*s0_slip/n_slip
!
!    *** Tolerances ***
 tol_in  = tol_inner*s0_slip
 tol_out = tol_outer*s0_slip 
!
 dgmaxc = 0
!    *** Error treatment ***
 iconv  = 0
 ising  = 0
!
!    *** Calculation of Fp_old(-1) ***
 invFp_old=Fp_old   
 call invert(invFp_old,3,0,0,det,3)
 if (det==0.0_pReal) then
    ising=1
    return
 endif
!
!    *** Calculation of A and T*0 (see Kalidindi) ***
! constitutive ÄÄÄ
 A=matmul(transpose(matmul(Fg_new,invFp_old)), matmul(Fg_new,invFp_old))
 call math_conv33to6((A-I3)/2,Estar0_v)
 Tstar0_v=matmul(Cslip_66,Estar0_v)
!
!    *** Calculation of Dslip (see Kalidindi) *** 
! constitutive ÄÄÄ
 do i=1,nslip
    mm=matmul(A,Sslip(i,:,:))
    mm1=(mm+transpose(mm))/2
    vv = math_33to6(mm1)
    Dslip(:,i)=matmul(Cslip_66,vv)
 enddo
!
!    *** Second level of iterative procedure: Resistences ***
 do iouter=1,nouter
!    *** First level of iterative procedure: Stresses ***
    do iinner=1,ninner
!
!    *** Calculation of gdot_slip ***    
! constitutive ÄÄÄ
        do i=1,nslip
            tau_slip(i)=dot_product(Tstar_v,Sslip_v(i,:))
        enddo
        call slip_rate(tau_slip,tauc_slip_new,gdot_slip,
                &     dgdot_dtaucslip)

!    *** Evaluation of Tstar and Gn (see Kalidindi) ***
        vv=0
        do i=1,nslip
            vv=vv-gdot_slip(i)*Dslip(:,i)
        enddo
        R1=Tstar_v-Tstar0_v-vv*dt  
        norm1=maxval(abs(R1))  
        if (norm1.LT.tol_in) goto 100
! 
!    *** Jacobi Calculation ***
        if (numerical==1) then
!    *** Perturbation method ***
        else
!    *** Analytical Calculation ***
            Jacobi=0
            do i=1,nslip
                do j=1,6
                    do k=1,6
                        Jacobi(j,k)=Jacobi(j,k)
                            &     +Dslip(j,i)*Sslip_v(i,k)*dgdot_dtaucslip(i)
                    enddo
                enddo
            enddo
            Jacobi=Jacobi*dt
            do i=1,6
                Jacobi(i,i)=1.0_pReal+Jacobi(i,i)
            enddo
        endif
!    *** End of the Jacobi calculation ***

!    *** Inversion of the Jacobi matrix ***
        invJacobi=Jacobi
        call invert(invJacobi,6,0,0,det,6)
        if (det==0.0_pReal) then
            do i=1,6
                Jacobi(i,i)=1.05d0*maxval(Jacobi(i,:))
            enddo
            invJacobi=Jacobi
            call invert(invJacobi,6,0,0,det,6)
            if (det==0.0_pReal) then
                ising=1
                return
            endif
        endif
        dTstar_v=matmul(invJacobi,R1)

!    *** Correction (see Kalidindi) ***
        do i=1,6
            if (abs(dTstar_v(i)).GT.crite) then
                dTstar_v(i)=sign(crite,dTstar_v(i))
            endif
        enddo
        Tstar_v=Tstar_v-dTstar_v
!
    enddo
    iconv=1
    return
!    *** End of the first level of iterative procedure ***

 100 continue

    call hardening(tauc_slip_new,gdot_slip,dtauc_slip)

!    *** Arrays of residuals ***
    R2=tauc_slip_new-tauc_slip_old-dtauc_slip*dt 
    norm2=maxval(abs(R2))
    if (norm2.LT.tol_out) goto 200
        tauc_slip_new=tauc_slip_old+dtauc_slip*dt
    enddo  
    iconv=2
    return
!    *** End of the second level of iterative procedure ***

 200 continue
!
 call plastic_vel_grad(dt,tau_slip,tauc_slip_new,Lp)
! 
!    *** Calculation of Fp(t+dt) (see Kalidindi) ***
 dLp=I3+Lp*dt
 Fp_new=matmul(dLp,Fp_old)
 call math_determ(Fp_new,det)
 Fp_new=Fp_new/det**(1.0_pReal/3.0_pReal)
!
!    *** Calculation of F*(t+dt) (see Kalidindi) ***
 invFp_new=Fp_new
 call invert(invFp_new,3,0,0,det,3)
 if (det==0.0_pReal) then
    ising=1
    return
 endif
 Fe=matmul(Fg_new,invFp_new)
!
!    *** Calculation of Estar ***
 Estar=0.5_pReal*(matmul(transpose(Fe),Fe)-I3)
 call CPFEM_conv33to6(Estar,Estar_v)
!
!    *** Calculation of the Cauchy stress ***
 call cauchy_stress(Estar_v,Fe,cs)
!
 return
 end
!
!
!
 end module