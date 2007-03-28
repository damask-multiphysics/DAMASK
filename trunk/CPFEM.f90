! last modified 28.03.07
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
 real(pReal), dimension (:,:,:),     allocatable :: CPFEM_stress_all
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_jacobi_all
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_ffn_all
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_ffn1_all
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_results
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_ini_ori
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_sigma_old
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_sigma_new
 real(pReal), dimension (:,:,:,:,:), allocatable :: CPFEM_Fp_old
 real(pReal), dimension (:,:,:,:,:), allocatable :: CPFEM_Fp_new
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_jaco_old
 integer(pInt) :: CPFEM_inc_old    = 0_pInt
 integer(pInt) :: CPFEM_subinc_old = 1_pInt
 integer(pInt) :: CPFEM_Nresults   = 3_pInt
 logical :: CPFEM_first_call = .true.

 CONTAINS

!***********************************************************************
!***    This routine checks for initialization, variables update and ***
!***    calls the actual material model                              ***
!***********************************************************************
 subroutine cpfem_general(ffn, ffn1, CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_dt, cp_en, CPFEM_in)
!
 use prec, only: pReal,pInt
! use CPFEM, only: CPFEM_ffn_all, CPFEM_ffn1_all, CPFEM_inc_old
! use IO, only: IO_init
 use constitutive, only: constitutive_state_old, constitutive_state_new
 implicit none
!
 real(pReal)   ffn(3,3), ffn1(3,3), CPFEM_dt
 integer(pInt) CPFEM_inc, CPFEM_subinc, CPFEM_cn, cp_en, CPFEM_in
!
! initialization step
 if (CPFEM_first_call) then
! three dimensional stress state ?
!    call IO_init()
    call mesh_init()
    call constitutive_init()
    call math_init()
    call CPFEM_init()
    CPFEM_first_call = .false.
 endif
! not a new increment
 if (CPFEM_inc==CPFEM_inc_old) then
! case of a new subincrement:update starting with subinc 2
     if (CPFEM_subinc > CPFEM_subinc_old) then
        CPFEM_sigma_old        = CPFEM_sigma_new
        CPFEM_Fp_old           = CPFEM_Fp_new
        constitutive_state_old = constitutive_state_new
        CPFEM_subinc_old       = CPFEM_subinc
    endif
! case of a new increment
 else
    CPFEM_sigma_old         = CPFEM_sigma_new
    CPFEM_Fp_old            = CPFEM_Fp_new
    constitutive_state_old  = constitutive_state_new
    CPFEM_inc_old           = CPFEM_inc
    CPFEM_subinc_old        = 1_pInt
 endif
!
! get cp element number for fe element number
 CPFEM_ffn_all(:,:,CPFEM_in, cp_en)  = ffn
 CPFEM_ffn1_all(:,:,CPFEM_in, cp_en) = ffn1
 call CPFEM_general_material(CPFEM_cn, CPFEM_dt, cp_en, CPFEM_in)
 return
 end subroutine


!***********************************************************************
!***    This routine allocates the arrays defined in module CPFEM    ***
!***    and initializes them                                         ***
!***********************************************************************
 subroutine CPFEM_init()
!
 use prec, only: pReal,pInt
! use math, only: math_I3
 use mesh
 use constitutive
!
 implicit none
!
 integer(pInt) i
!
!    *** mpie.marc parameters ***
 allocate(CPFEM_ffn_all(3,3,mesh_maxNips,mesh_NcpElems))
 allocate(CPFEM_ffn1_all(3,3,mesh_maxNips,mesh_NcpElems))
 allocate(CPFEM_stress_all(6,mesh_maxNips,mesh_NcpElems))
 allocate(CPFEM_jacobi_all(6,6,mesh_maxNips,mesh_NcpElems))
 CPFEM_ffn_all    = 0.0_pReal
 CPFEM_ffn1_all   = 0.0_pReal
 CPFEM_stress_all = 0.0_pReal
 CPFEM_jacobi_all = 0.0_pReal
!
!    *** User defined results !!! MISSING incorporate consti_Nresults ***
 allocate(CPFEM_results(CPFEM_Nresults+constitutive_maxNresults,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_results = 0.0_pReal
!
!    *** Second Piola-Kirchoff stress tensor at (t=t0) and (t=t1) ***
 allocate(CPFEM_sigma_old(6,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(CPFEM_sigma_new(6,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_sigma_old = 0.0_pReal
 CPFEM_sigma_new = 0.0_pReal
!
!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***  
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_Fp_old = 0.0_pReal 
 CPFEM_Fp_new = 0.0_pReal 
!    
!    *** Old jacobian (consistent tangent) ***
 allocate(CPFEM_jaco_old(6,6,mesh_maxNips,mesh_NcpElems))
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
 write(6,*) 'CPFEM_sigma_old:     ', shape(CPFEM_sigma_old)
 write(6,*) 'CPFEM_sigma_new:     ', shape(CPFEM_sigma_new)
 write(6,*) 'CPFEM_Fp_old:        ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:        ', shape(CPFEM_Fp_new)
 write(6,*) 'CPFEM_jaco_old:      ', shape(CPFEM_jaco_old)
 write(6,*)
 call flush(6)
 return
 end subroutine  
!
!
 subroutine CPFEM_general_material(&
     CPFEM_cn,&       ! Cycle number
     CPFEM_dt,&       ! Time increment (dt)
     cp_en,&          ! Element number
     CPFEM_in)        ! Integration point number
!***********************************************************************
!***    This routine calculates the material behaviour    ***
!***********************************************************************
 use prec, only: pReal,pInt
! use IO, only: IO_error
 use math
 use mesh
 use constitutive
!
 implicit none
!
!    *** Definition of variables ***
!    *** Subroutine parameters ***
 real(pReal)   CPFEM_dt
 integer(pInt) CPFEM_cn, cp_en ,CPFEM_in
!    *** Local variables ***
 real(pReal)   vf, cs(6), cd(6,6), CPFEM_d(6,6), CPFEM_s(6)
 integer(pInt) jpara,nori, iori, ising, icut, iconv, CPFEM_en
!    *** Numerical parameters ***
!    *** How often the jacobian is recalculated ***
 integer (pInt), parameter :: ijaco = 5_pInt
!    *** Reference shear rate for the calculation of CPFEM_timefactor ***
 real (pReal), parameter :: dgs = 0.01_pReal
!
!    *** Flag for recalculation of jacobian ***
 jpara = 1_pInt      
! get number of grains from cp element number and integration point number
 nori = constitutive_Ngrains(CPFEM_in,cp_en)  !ÄÄÄ
! 
 CPFEM_en = mesh_element(1,cp_en) ! remap back to FE id
!
 CPFEM_s=0
 CPFEM_d=0
!
!    *** Loop over all the components ***
 do iori=1,nori
!
!    *** Initialization of the matrices for t=t0 ***
! data from constitutive?
    vf = constitutive_matVolFrac(iori,CPFEM_in,cp_en)*constitutive_texVolFrac(iori,CPFEM_in,cp_en) !ÄÄÄ

!    *** Calculation of the solution at t=t1 ***
!	QUESTION use the mod() as flag parameter in the call ??
    if (mod(CPFEM_cn,ijaco)==0) then !ÄÄÄ
        call CPFEM_stress(cs, cd, CPFEM_dt,cp_en,CPFEM_in, iori, ising, icut, iconv, 1_pInt)
!    ***   Evaluation of ising      ***
!    *** ising=2 => singular matrix in jacobi calculation ***
!    ***      => use old jacobi        ***
        if (ising==2) jpara=0
!    *** Calculation of the consistent tangent ***
        CPFEM_d=CPFEM_d+vf*cd
    else
       call CPFEM_stress(cs, cd, CPFEM_dt,cp_en,CPFEM_in, iori, ising, icut, iconv, 0_pInt)
       jpara=0
    endif
!    *** Cases of unsuccessful calculations *** 
!    ***   Evaluation of ising  ***
!    *** ising!=0 => singular matrix ***
    if (ising==1) then
        write(6,*) 'Singular matrix!'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(700)
!        CPFEM_timefactor=1.e5_pReal
        return
    endif
!    ***   Evaluation of icut   ***
!    *** icut!=0 => too many cutbacks ***
    if (icut==1) then
        write(6,*) 'Too many cutbacks'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(600)
!        CPFEM_timefactor=1.e5_pReal
        return
    endif
!    ***  Evaluation of iconv ***
!    *** iconv!=0 => no convergence ***
    if (iconv==1) then
        write(6,*) 'Inner loop did not converge!'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(600)
!       CPFEM_timefactor=1.e5_pReal
        return
    else if (iconv==2) then
        write(6,*) 'Outer loop did not converge!'
        write(6,*) 'Integration point: ',CPFEM_in
        write(6,*) 'Element:           ',CPFEM_en
        call IO_error(600)
!        CPFEM_timefactor=1.e5_pReal
        return
    endif
!    *** Evaluation of the average Cauchy stress ***    
    CPFEM_s=CPFEM_s+vf*cs
 enddo
!    *** End of the loop over the components ***

!    *************************************
!    *** End of the CP-FEM Calculation *** 
!    *************************************
!    *** Store the new stress ***     
 CPFEM_stress_all(:,CPFEM_in,cp_en)=CPFEM_s
!    *** Store the new jacobian ***     
 if (jpara/=0) CPFEM_jaco_old(:,:,CPFEM_in,cp_en)=CPFEM_d
 return
 end subroutine
!
!
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
     isjaco)            ! flag whether to calculate Jacoby matrix
!********************************************************************
! This routine calculates the stress for a single component
! and manages the independent time incrmentation
!********************************************************************
 use prec, only: pReal,pInt
 use constitutive, only: constitutive_Nstatevars, constitutive_state_old, constitutive_state_new, constitutive_Nresults,&
                         constitutive_results
 implicit none
!
!    *** Definition of variables ***
!    *** Subroutine parameters ***
 real(pReal)   cs(6), cd(6,6), CPFEM_dt
 integer(pInt) cp_en ,CPFEM_in, iori, ising, icut, iconv, isjaco
!    *** Local variables ***
 real(pReal)   Fp_old(3,3), Fp_new(3,3), state_old(constitutive_Nstatevars(iori, CPFEM_in, cp_en))
 real(pReal)   state_new(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), Tstar_v(6), CPFEM_ffn(3,3), CPFEM_ffn1(3,3)
 real(pReal)   Tstar_v_h(6), state_new_h(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), phi1, PHI, phi2, dt_i
 real(pReal)   delta_Fg(3,3), Fg_i(3,3), state_new_i(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), time
 integer(pInt) jcut
!    *** Numerical parameters ***
 integer(pInt), parameter :: ncut=7_pInt
!
 icut=0
!
!    *** Initialization of the matrices for t=t0 ***
    Fp_old     = CPFEM_Fp_old(:,:,iori,CPFEM_in,cp_en)
    Fp_new     = 0.0_pReal
    state_old  = constitutive_state_old(:,iori,CPFEM_in,cp_en)
    state_new  = state_old
    Tstar_v    = CPFEM_sigma_old(:,iori,CPFEM_in,cp_en)
    CPFEM_ffn  = CPFEM_ffn_all(:,:,CPFEM_in,cp_en)
    CPFEM_ffn1 = CPFEM_ffn1_all(:,:,CPFEM_in,cp_en)
!
!    *** First attempt to calculate Tstar and tauc with initial timestep ***
! save copies of Tstar_v and state_new
 Tstar_v_h = Tstar_v
 state_new_h = state_new
 call CPFEM_stress_int(cs, cd, CPFEM_dt, cp_en,CPFEM_in, iori, ising, icut, iconv, isjaco, phi1, PHI, phi2,&
                       CPFEM_ffn, CPFEM_ffn1,Fp_old,Fp_new,state_old, state_new, Tstar_v)
 if ((iconv==0).AND.(ising==0)) then
!    *** Update the differents matrices for t=t1 ***
    CPFEM_Fp_new(:,:,iori,CPFEM_in,cp_en)          = Fp_new
    constitutive_state_new(:,iori,CPFEM_in,cp_en) = state_new
    CPFEM_sigma_new(:,iori,CPFEM_in,cp_en)         = Tstar_v
!    *** Update the results plotted in MENTAT ***
    CPFEM_results(1,iori,CPFEM_in,cp_en) = phi1
    CPFEM_results(2,iori,CPFEM_in,cp_en) = PHI
    CPFEM_results(3,iori,CPFEM_in,cp_en) = phi2
    CPFEM_results(4:3+constitutive_Nresults(iori,CPFEM_in,cp_en),iori,CPFEM_in,cp_en)=&
            constitutive_results(1:constitutive_Nresults(iori,CPFEM_in,cp_en),iori,CPFEM_in,cp_en)!ÄÄÄÄ
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
    call CPFEM_stress_int(cs, cd, time, cp_en,CPFEM_in, iori, ising, icut, iconv, isjaco, phi1, PHI, phi2,&
                       CPFEM_ffn, Fg_i,Fp_old,Fp_new,state_old, state_new_i, Tstar_v)
    if ((iconv==0).AND.(ising==0)) then
        time=time+dt_i
        Fg_i=Fg_i+delta_Fg
        Tstar_v_h=Tstar_v
        state_new_h=state_new_i
    else
        jcut=jcut+1_pInt
        if (jcut>ncut) then
            icut=1_pInt
            return
        endif
        dt_i=0.5_pReal*dt_i
        time=time-dt_i
        delta_Fg=0.5_pReal*delta_Fg
        Fg_i=Fg_i-delta_Fg
        Tstar_v=Tstar_v_h
        state_new_i=state_new_h
    endif
 enddo
!
!    *** Final calculation of stress and resistences with full timestep ***
 state_new=state_new_i
 call CPFEM_stress_int(cs, cd, CPFEM_dt, cp_en,CPFEM_in, iori, ising, icut, iconv, isjaco, phi1, PHI, phi2,&
                       CPFEM_ffn, CPFEM_ffn1,Fp_old,Fp_new,state_old, state_new, Tstar_v)
!    *** Update the differents matrices for t=t1 ***
 CPFEM_Fp_new(:,:,iori,CPFEM_in,cp_en)          = Fp_new
 constitutive_state_new(:,iori,CPFEM_in,cp_en) = state_new
 CPFEM_sigma_new(:,iori,CPFEM_in,cp_en)         = Tstar_v
!    *** Update the results plotted in MENTAT ***
 CPFEM_results(1,iori,CPFEM_in,cp_en) = phi1
 CPFEM_results(2,iori,CPFEM_in,cp_en) = PHI
 CPFEM_results(3,iori,CPFEM_in,cp_en) = phi2
 return
 end subroutine
!
!
 subroutine CPFEM_stress_int(&
     cs,&         ! Cauchy stress vector
     dcs_de,&     ! Consistent tangent
     dt,&         ! Time increment
     cp_en,&      ! Element number
     CPFEM_in,&   ! Integration point number
     iori,&       ! number of orintation
     ising,&      ! flag for singular matrix
     icut,&       ! flag for too many cut backs
     iconv,&      ! flag for non convergence
     isjaco,&     ! flag whether to calculate Jacoby matrix
     phi1,&       ! Euler angle
     PHI,&        ! Euler angle
     phi2,&       ! Euler angle
     Fg_old,&     ! Old global deformation gradient
     Fg_new,&     ! New global deformation gradient
     Fp_old,&     ! Old plastic deformation gradient
     Fp_new,&     ! New plastic deformation gradient
     state_old,&  ! Old state variable array
     state_new,&  ! New state variable array
     Tstar_v)     ! Second Piola-Kirschoff stress tensor
!********************************************************************
! This routine calculates the stress for a single component
! it is based on the paper by Kalidindi et al.:
! J. Mech. Phys, Solids Vol. 40, No. 3, pp. 537-569, 1992
! it is modified to use anisotropic elasticity matrix
!********************************************************************
 use prec, only: pReal,pInt
 use constitutive, only: constitutive_Nstatevars
 use math, only: math_Mandel6to33
 implicit none
!
!    *** Definition of variables ***
!    *** Subroutine parameters ***
 integer(pInt) cp_en, CPFEM_in, iori, ising, icut, iconv, isjaco
 real(pReal)   cs(6), dcs_de(6,6), dt, phi1, PHI, phi2, Fg_old(3,3), Fg_new(3,3)
 real(pReal)   Fp_old(3,3), Fp_new(3,3), state_old(constitutive_Nstatevars(iori, CPFEM_in, cp_en))
 real(pReal)   state_new(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), Tstar_v(6)
!    *** Local variables ***
 integer(pInt) ic     
 real(pReal) Fe(3,3), R(3,3), U(3,3), dev(6), dF(3,3), Fg2(3,3), sgm2(6)
 real(pReal) state2(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), Fp2(3,3), cs1(6)
!    *** Numerical parameters ***
 real(pReal), parameter :: pert_ct=1.0e-5_pReal  
!    *** Error treatment ***
 iconv  = 0
 ising  = 0   

!    *********************************************
!    ***  Calculation of the new Cauchy stress ***
!    *********************************************

!    *** Call Newton-Raphson method ***
 call NEWTON_RAPHSON(dt,cp_en,CPFEM_in,iori,Fg_old,Fg_new,Fp_old,Fp_new,Fe,state_old,state_new,Tstar_v,cs,iconv,ising)
! 
!    *** Calculation of the new orientation ***
 call math_pDecomposition(Fe,U,R,ising)
 if (ising==1) then
    return
 endif
 call math_RtoEuler(transpose(R),phi1,PHI,phi2)
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
    dF=matmul(math_Mandel6to33(dev),Fg_old) 
    Fg2=Fg_new+dF
    sgm2=Tstar_v
    state2=state_new

!    *** Calculation of the perturbated Cauchy stress ***
    call NEWTON_RAPHSON(dt,cp_en,CPFEM_in,iori,Fg_old,Fg2,Fp_old,Fp2,Fe,state_old,state2,sgm2,cs1,iconv,ising)
!  
!    *** Consistent tangent ***
    dcs_de(:,ic)=(cs1-cs)/pert_ct     
 enddo
!
 return
 end subroutine
!
!
 subroutine NEWTON_RAPHSON(&
     dt,&
     cp_en,&            ! Element number
     CPFEM_in,&         ! Integration point number
     iori,&             ! number of orintation
     Fg_old,&
     Fg_new,&
     Fp_old,&
     Fp_new,&
     Fe,&
     state_old,&
     state_new,&
     Tstar_v,&
     cs,&
     iconv,&
     ising)
!***********************************************************************
!***        NEWTON-RAPHSON Calculation                               ***
!***********************************************************************
 use prec, only: pReal,pInt
 use constitutive, only: constitutive_Nstatevars, constitutive_HomogenizedC, constitutive_dotState
 use math
 implicit none
!    *** Definition of variables ***
!    *** Subroutine parameters ***
 integer(pInt) cp_en, CPFEM_in, iori, iconv, ising
 real(pReal) dt,Fg_old(3,3),Fg_new(3,3),Fp_old(3,3),Fp_new(3,3), Fe(3,3)
 real(pReal) state_old(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), state_new(constitutive_Nstatevars(iori, CPFEM_in, cp_en)) 
 real(pReal) Tstar_v(6), cs(6)
!    *** Local variables ***
 real(pReal) tol_in, tol_out, invFp_old(3,3), det, A(3,3), C_66(6,6), Lp(3,3), dLp(3,3,3,3)
 real(pReal) I3tLp(3,3), help(3,3), help1(3,3,3,3), Tstar0_v(6), R1(6), R1s(6), norm1, tdLp(3,3)
 real(pReal) dstate(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), R2(constitutive_Nstatevars(iori, CPFEM_in, cp_en))
 real(pReal) R2s(constitutive_Nstatevars(iori, CPFEM_in, cp_en)), norm2, invFp_new(3,3), Estar(3,3)
 real(pReal) Estar_v(6), Jacobi(6,6), invJacobi(6,6), dTstar_v(6), help2(6,6)
 integer(pInt) iouter, iinner , dummy, err, i, j, k, l, m
!    *** Numerical parameters ***
 integer(pInt), parameter :: nouter    = 50_pInt
 real(pReal),   parameter :: tol_outer = 1.0e-4_pReal
 integer(pInt), parameter :: ninner    = 2000_pInt
 real(pReal),   parameter :: tol_inner = 1.0e-3_pReal   
 real(pReal),   parameter :: crite     = 1.0e-1_pReal 

! crite=eta*constitutive_s0_slip/constitutive_n_slip !ÄÄÄ
!
!    *** Tolerances ***
! tol_in  = tol_inner*s0_slip !ÄÄÄ
! tol_out = tol_outer*s0_slip !ÄÄÄ
!
!    *** Error treatment ***
 iconv  = 0
 ising  = 0
!
! initialize new state
 state_new=state_old
!    *** Calculation of Fp_old(-1) ***
 call invert3x3(Fp_old, invFp_old, det, err) !ÄÄÄ
 if (err==1_pInt) then
    ising=1
    return
 endif
!
!    *** Calculation of A and T*0 (see Kalidindi) ***
 A = matmul(Fg_new,invFp_old)  ! actually Fe
 A = matmul(transpose(A), A)
 C_66=constitutive_HomogenizedC(iori, CPFEM_in, cp_en) !ÄÄÄ
 Tstar_v=matmul(C_66, math_Mandel33to6(A-math_I3))  ! fully elastic guess
! QUESTION follow former plastic slope to guess better?
!
!    *** Second level of iterative procedure: Resistences ***
 do iouter=1,nouter
!    *** First level of iterative procedure: Stresses ***
    do iinner=1,ninner
!
!    *** Calculation of gdot_slip ***    
        call constitutive_LpAndItsTangent(Tstar_v, iori, CPFEM_in, cp_en, Lp, dLp)
        I3tLp  = math_I3-dt*Lp
        help=matmul(transpose(I3tLp),matmul(A, I3tLp))-math_I3
        Tstar0_v = 0.5_pReal * matmul(C_66, math_Mandel33to6(help))
        R1=Tstar_v-Tstar0_v
        if (maxval(abs(R1/maxval(abs(Tstar_v)))) < tol_inner) goto 100
! 
!    *** Jacobi Calculation: dRes/dTstar ***
        help=matmul(A, I3tLp)
        help1=0
        do i=1,3
          do j=1,3
            do k=1,3
              do l=1,3
                do m=1,3
                  help1(i,j,k,l)=help1(i,j,k,l)+help(i,m)*dLp(m,j,k,l)+help(j,m)*dLp(m,i,l,k)
                enddo
              enddo
            enddo
          enddo
        enddo
        help2=math_Mandel3333to66(help1)
        Jacobi= 0.5_pReal*matmul(C_66, help2) + math_identity2nd(6)
        call math_invert6x6(Jacobi, invJacobi, dummy, err) !ÄÄÄ
        if (err==1_pInt) then
            forall (i=1:6) Jacobi(i,i)=1.05d0*maxval(Jacobi(i,:)) ! regularization
            call math_invert6x6(Jacobi, invJacobi, dummy, err)
            if (err==1_pInt) then  ! sorry, can't help here!!
                ising=1
                return
            endif
        endif
        dTstar_v=matmul(invJacobi,R1)  ! correction to Tstar

!    *** Correction (see Kalidindi) ***
        forall(i=1:6, abs(dTstar_v(i)) > crite*maxval(abs(Tstar_v))) &
          dTstar_v(i) = sign(crite*maxval(abs(Tstar_v)),dTstar_v(i))

        Tstar_v=Tstar_v-dTstar_v
!
    enddo
    iconv=1
    return
!    *** End of the first level of iterative procedure ***

100 dstate=dt*constitutive_dotState(Tstar_v, iori, CPFEM_in, cp_en)
!    *** Arrays of residuals ***
    R2=state_new-state_old-dstate
    R2s=0.0_pReal
    forall(i=1:constitutive_Nstatevars(iori, CPFEM_in, cp_en), state_new(i)/=0.0_pReal) R2s(i)=R2(i)/state_new(i)
    if (maxval(dabs(R2s)) < tol_outer) goto 200
    state_new=state_old+dstate
 enddo  
 iconv=2
 return
!    *** End of the second level of iterative procedure ***

!    *** Calculation of Fp(t+dt) (see Kalidindi) ***
200 invFp_new=matmul(Fp_old, I3tLp)
 call math_invert3x3(invFp_new, Fp_new, det, err) !ÄÄÄ
 if (err==1_pInt) then
    ising=1
    return
 endif
 Fp_new=Fp_new/math_det3x3(Fp_new)**(1.0_pReal/3.0_pReal)
!
!    *** Calculation of F*(t+dt) (see Kalidindi) ***
 Fe=matmul(Fg_new,invFp_new)
!
!    *** Calculation of the Cauchy stress ***
!	QUESTION seems to need Tstar, not Estar..??
 cs = CPFEM_cauchy_stress(Tstar_v,Fe)
!
 return
 end subroutine
!
 function CPFEM_cauchy_stress(PK_v, Fe)
!***********************************************************************
!***        Cauchy stress calculation                               ***
!***********************************************************************
 use prec, only: pReal,pInt
 use math, only: math_Mandel33to6,math_Mandel6to33,math_det3x3
 implicit none
!    *** Subroutine parameters ***
 real(pReal) PK_v(6), Fe(3,3), CPFEM_cauchy_stress(6)

 CPFEM_cauchy_stress = math_Mandel33to6(matmul(matmul(Fe,math_Mandel6to33(PK_v)),transpose(Fe))/math_det3x3(Fe))
 end function
 end module