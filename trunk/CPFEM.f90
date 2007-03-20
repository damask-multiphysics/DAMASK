
!    ---------------------------
 MODULE CPFEM     
!    ---------------------------
!    *** CPFEM engine ***

 use prec, only: pRe,pIn
 implicit none

!    ****************************************************************
!    *** General variables for the material behaviour calculation ***
!    ****************************************************************  
 real(pRe), allocatable :: CPFEM_stress_all   (:,:,:)
 real(pRe), allocatable :: CPFEM_jacobi_all   (:,:,:,:)
 real(pRe), allocatable :: CPFEM_results    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_thickness    (:,:)
 real(pRe), allocatable :: CPFEM_ini_ori    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_sigma_old    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_sigma_new    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_Fp_old    (:,:,:,:,:)
 real(pRe), allocatable :: CPFEM_Fp_new    (:,:,:,:,:)
 real(pRe), allocatable :: CPFEM_tauc_slip_old(:,:,:,:)
 real(pRe), allocatable :: CPFEM_tauc_slip_new(:,:,:,:)
 real(pRe), allocatable :: CPFEM_g_old    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_g_new    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_jaco_old    (:,:,:,:)
 real(pRe), allocatable :: CPFEM_mat     (:,:)

 CONTAINS

!***********************************************************************
!***    This routine allocates the arrays defined in module mpie    ***
!***    and initializes them          ***
!***********************************************************************
 subroutine ALLOCATION(mpie_numel,mpie_nip)

 use prec, only: pRe,pIn
 use IO, only: _error
 use math
 use mesh
 use constitutive
 
 implicit none

 integer(pIn) i

!    *** mpie.marc parameters ***
 allocate(CPFEM_stress_all(6,mesh_Nelems,mesh_Nips))
 allocate(CPFEM_jacobi_all(6,6,mesh_Nelems,mesh_Nips))
 CPFEM_stress_all=0.0_pRe
 CPFEM_jacobi_all=0.0_pRe

!    *** User defined results ***
 allocate(CPFEM_results(constitutive_Nresults,
     &        constitutive_maxNgrains,
     &        mesh_Nelems,mesh_Nips))
 CPFEM_results=0.0_pRe

!    *** Relative sheet thickness ***
 allocate(CPFEM_thickness(mesh_Nelems,mesh_Nips))
 CPFEM_thickness=0.0_pRe

!    *** Initial orientations ***
 allocate(CPFEM_ini_ori(3,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))
 CPFEM_ini_ori=0.0_pRe

!    *** Second Piola-Kirchoff stress tensor at (t=t0) and (t=t1) ***
 allocate(CPFEM_sigma_old(6,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))
 allocate(CPFEM_sigma_new(6,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))
 CPFEM_sigma_old=0.0_pRe
 CPFEM_sigma_new=0.0_pRe

!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***  
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))
 CPFEM_Fp_old=0.0_pRe 
 CPFEM_Fp_new=0.0_pRe 
 do i=1,3
    CPFEM_Fp_old(i,i,:,:,:)=1.0_pRe 
    CPFEM_Fp_new(i,i,:,:,:)=1.0_pRe
 enddo 
     
! QUESTION: would it be wise to outsource these to _constitutive_ ??
!    *** Slip resistances at (t=t0) and (t=t1) ***  
 allocate(CPFEM_tauc_slip_old(nslip,constitutive_maxNgrains,mesh_Nelems,
     &        mesh_Nips))   
 allocate(CPFEM_tauc_slip_new(nslip,constitutive_maxNgrains,mesh_Nelems,
     &        mesh_Nips))
 CPFEM_tauc_slip_old=0.0_pRe
 CPFEM_tauc_slip_new=0.0_pRe

!    *** Cumulative shear at (t=t0) and (t=t1) *** 
! QUESTION which nslip to use here ?!?
 allocate(CPFEM_g_old(nslip,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))  
 allocate(CPFEM_g_new(nslip,constitutive_maxNgrains,mesh_Nelems,mesh_Nips))
 CPFEM_g_old=0.0_pRe
 CPFEM_g_new=0.0_pRe

!    *** Old jacobian (consistent tangent) ***
 allocate(CPFEM_jaco_old(6,6,mesh_Nelems,mesh_Nips))

!    *** Output to MARC output file ***
 write(6,*)
 write(6,*) 'Arrays allocated:'
 write(6,*) 'CPFEM_stress_all:    ', shape(CPFEM_stress_all)
 write(6,*) 'CPFEM_jacobi_all:    ', shape(CPFEM_jacobi_all)
 write(6,*) 'CPFEM_results:    ', shape(CPFEM_results)
 write(6,*) 'CPFEM_thickness:    ', shape(CPFEM_thickness)
 write(6,*) 'CPFEM_ini_ori:    ', shape(CPFEM_ini_ori)
 write(6,*) 'CPFEM_sigma_old:    ', shape(CPFEM_sigma_old)
 write(6,*) 'CPFEM_sigma_new:    ', shape(CPFEM_sigma_new)
 write(6,*) 'CPFEM_Fp_old:    ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:    ', shape(CPFEM_Fp_new)
 write(6,*) 'CPFEM_tauc_slip_old: ', shape(CPFEM_tauc_slip_old)
 write(6,*) 'CPFEM_tauc_slip_new: ', shape(CPFEM_tauc_slip_new)
 write(6,*) 'CPFEM_g_old:    ', shape(CPFEM_g_old)
 write(6,*) 'CPFEM_g_new:    ', shape(CPFEM_g_new)
 write(6,*) 'CPFEM_jaco_old:    ', shape(CPFEM_jaco_old)
 write(6,*)
 call flush(6)
 return
 end     
 

 subroutine CPFEM_general_material(
     & CPFEM_s,   ! Stress vector
     & CPFEM_d,   ! Jacobi matrix (consistent tangent)
     & CPFEM_ndi,   ! Dimension
     & CPFEM_ffn,   ! Deformation gradient at begin of increment
     & CPFEM_ffn1,  ! Deformation gradient at end of increment
     & CPFEM_inc,   ! Increment number
     & CPFEM_subinc,  ! Subincrement number
     & CPFEM_cn,   ! Cycle number
     & CPFEM_tinc,  ! Time increment (dt)
     & CPFEM_timefactor, ! Factor for timestep correction
! & mesh_Nelems,   ! Number of elements in mesh
! & CPFEM_nip,   ! Maximum number of integration points per element
     & CPFEM_en,   ! Element number
     & CPFEM_in,   ! Integration point number
     & CPFEM_mn,   ! Material number
     & CPFEM_dimStress  ! Dimension of stress/strain vector
     &)
!***********************************************************************
!***    This routine calculates the material behaviour    ***
!***********************************************************************
 use prec, only: pRe,pIn
 use IO, only _error
 use math
 use mesh
 use constitutive
 
 implicit none

!    *** Definition of variables ***
 integer(pIn) CPFEM_ndi,CPFEM_inc,CPFEM_subinc,CPFEM_cn,
     &   CPFEM_en,CPFEM_in,CPFEM_mn,CPFEM_dimStress
 real(pRe) CPFEM_timefactor,CPFEM_tinc,CPFEM_s(CPFEM_dimStress),
     &     CPFEM_d(CPFEM_dimStress,CPFEM_dimStress),
     &     CPFEM_ffn(3,3),CPFEM_ffn1(3,3)
! QUESTION which nslip to use?
 real(pRe) Fp_old(3,3),tauc_slip_old(nslip),
     &     tauc_slip_new(nslip),g_old(nslip),
     &     g_new(nslip),Tstar_v(6),
     &     Fp_new(3,3),cs(6),phi1mis(2),PHImis(2),phi2mis(2),
     &     cd(6,6),ori_mat(3,3),hh6(6,6)
 integer(pIn) jpara,nori
 real(pRe) phi1,PHI,phi2,scatter,vf,alpha1,alpha2,beta1,
     &      beta2,phi1_s,PHI_s,phi2_s,p10,P0,p20,p11,P1,p21,
     &      dgmax,dgmaxc,orimis
 integer(pIn) i,iori,iconv,ising,icut
!    *** Numerical parameters ***
!    *** How often the jacobian is recalculated ***
 integer (pIn), parameter :: ijaco=1_pIn
!    *** Reference shear rate for the calculation of CPFEM_timefactor ***
 real (pRe), parameter :: dgs=0.01_pRe

!    *** Initialization step ***   
 if (CPFEM_first_call==1_pIn) then
    call INITIALIZATION(mesh_Nelems,CPFEM_nip)
    CPFEM_first_call=0_pIn
 endif
!    *** Case of a new increment *** 
 if (CPFEM_inc.NE.CPFEM_inc_old) then
    CPFEM_sigma_old=CPFEM_sigma_new
    CPFEM_Fp_old=CPFEM_Fp_new
    CPFEM_tauc_slip_old=CPFEM_tauc_slip_new
    CPFEM_g_old=CPFEM_g_new
    CPFEM_inc_old=CPFEM_inc
    CPFEM_subinc_old=1_pIn
    CPFEM_timefactor_max=0.0_pRe
 endif
!    *** case of a new subincrement:update starting with subinc 2 ***
 if (CPFEM_subinc.GT.CPFEM_subinc_old) then
    CPFEM_sigma_old=CPFEM_sigma_new
    CPFEM_Fp_old=CPFEM_Fp_new
    CPFEM_tauc_slip_old=CPFEM_tauc_slip_new
    CPFEM_g_old=CPFEM_g_new
    CPFEM_subinc_old=CPFEM_subinc
 endif
!    *** Flag for recalculation of jacobian ***
 jpara=1_pIn      

!    ************************************
!    ***  Orientation initialization  *** 
!    ************************************
!    *** Number of components per state ***  
 nori=CPFEM_mat(CPFEM_mn,1)
 if (CPFEM_inc==0_pIn) then
!    *** Three dimensional stress state ***
    if (CPFEM_ndi.NE.3_pIn) then
  call CPFEM_error(300)
    endif

    if ((CPFEM_en==1_pIn).AND.(CPFEM_in==1_pIn)) then
  write(6,*) 'MPIE Material Routine Ver. 0.1 by L. Hantcherli'
  write(6,*)
  write(6,*) 'Orientation initialization'
  call flush(6)
    endif

    i=1
    do while (i.LE.nori)
!    *** Direct ODF sampling ***
  if (CPFEM_mat(CPFEM_mn,2)==2) then
     call CPFEM_odf_ori(CPFEM_cko(CPFEM_mn,:,:,:,:), 
     &       CPFEM_odfmax(CPFEM_mn),phi1,PHI,phi2)
  else
!    *** Gauss/Spherical component ***
     if (CPFEM_mat(CPFEM_mn,7*i-4)==1) then
   phi1=CPFEM_mat(CPFEM_mn,7*i-3)
   PHI=CPFEM_mat(CPFEM_mn,7*i-2)
   phi2=CPFEM_mat(CPFEM_mn,7*i-1)
   scatter=CPFEM_mat(CPFEM_mn,7*i+1)
!    *** Random orientation to this component to represent ***
!    ***  random fraction of texture using halton series   ***
   if (phi1==400.0) then
      call CPFEM_halton_ori(phi1,PHI,phi2,scatter)
!    *** ELSE modify orientation to represent gauss distribution ***
   else if (scatter.GT.0.1) then
      call CPFEM_gauss(phi1,PHI,phi2,scatter)
   endif
!    *** Fiber component ***    
     else if (CPFEM_mat(CPFEM_mn,7*i-4)==2) then
   alpha1=CPFEM_mat(CPFEM_mn,7*i-3)
   alpha2=CPFEM_mat(CPFEM_mn,7*i-2)
   beta1=CPFEM_mat(CPFEM_mn,7*i-1)
   beta2=CPFEM_mat(CPFEM_mn,7*i)
   scatter=CPFEM_mat(CPFEM_mn,7*i+1)
!    *** Random orientation to this component to represent ***
!    ***  random fraction of texture using random numbers  ***
   if (alpha1==400.0) then
      call CPFEM_random_ori(phi1,PHI,phi2,scatter)
!    *** ELSE calculate orientation to represent fiber component ***
   else if (scatter.GT.0.1) then
      call CPFEM_fiber(alpha1,alpha2,beta1,beta2,
     &         scatter,phi1,PHI,phi2)
   endif
     else
   call CPFEM_error(510)
     endif
  endif
  CPFEM_ini_ori(1,i,CPFEM_en,CPFEM_in)=phi1
  CPFEM_ini_ori(2,i,CPFEM_en,CPFEM_in)=PHI
  CPFEM_ini_ori(3,i,CPFEM_en,CPFEM_in)=phi2
!    *** Orientation matrix ***    
  call CPFEM_euldreh(phi1,PHI,phi2,ori_mat)
  CPFEM_Fp_old(:,:,i,CPFEM_en,CPFEM_in)=ori_mat
  i=i+1

!    *** If symmetric component, creation of additional three orientations ***
  if (CPFEM_mat(CPFEM_mn,2)==1) then
!    *** First one ***
     phi1_s=180.0_pRe-phi1
     if (phi1_s.LT.0.0_pRe) phi1_s=phi1_s+360.0_pRe
     PHI_s=180.0_pRe-PHI
     if (PHI_s.LT.0.0_pRe) PHI_s=PHI_s+360.0_pRe
     phi2_s=phi2+180.0_pRe
     if (phi2_s.GT.360.0_pRe) phi2_s=phi2_s-360.0_pRe
     CPFEM_ini_ori(1,i,CPFEM_en,CPFEM_in)=phi1_s
     CPFEM_ini_ori(2,i,CPFEM_en,CPFEM_in)=PHI_s
     CPFEM_ini_ori(3,i,CPFEM_en,CPFEM_in)=phi2_s
!    *** Orientation matrix for initial orientation ***
     call CPFEM_euldreh(phi1_s,PHI_s,phi2_s,ori_mat)
     CPFEM_Fp_old(:,:,i,CPFEM_en,CPFEM_in)=ori_mat
     i=i+1
!    *** Second one ***
     phi1_s=360.0_pRe-phi1
     PHI_s=180.0_pRe-PHI
     if (PHI_s.LT.0.0_pRe) PHI_s=PHI_s+360.0_pRe
     phi2_s=phi2+180.0_pRe
     if (phi2_s.GT.360.0_pRe) phi2_s=phi2_s-360.0_pRe
     CPFEM_ini_ori(1,i,CPFEM_en,CPFEM_in)=phi1_s
     CPFEM_ini_ori(2,i,CPFEM_en,CPFEM_in)=PHI_s
     CPFEM_ini_ori(3,i,CPFEM_en,CPFEM_in)=phi2_s
!    *** Orientation matrix for initial orientation ***
     call CPFEM_euldreh(phi1_s,PHI_s,phi2_s,ori_mat)
     CPFEM_Fp_old(:,:,i,CPFEM_en,CPFEM_in)=ori_mat
     i=i+1
!    *** Third one ***
     phi1_s=phi1+180.0_pRe
     if (phi1_s.GT.360.0_pRe) phi1_s=phi1_s-360.0_pRe
     PHI_s=PHI
     phi2_s=phi2
     CPFEM_ini_ori(1,i,CPFEM_en,CPFEM_in)=phi1_s
     CPFEM_ini_ori(2,i,CPFEM_en,CPFEM_in)=PHI_s
     CPFEM_ini_ori(3,i,CPFEM_en,CPFEM_in)=phi2_s
!    *** Orientation matrix for initial orientation ***
     call CPFEM_euldreh(phi1_s,PHI_s,phi2_s,ori_mat)
     CPFEM_Fp_old(:,:,i,CPFEM_en,CPFEM_in)=ori_mat
     i=i+1
  else if ((CPFEM_mat(CPFEM_mn,2).NE.0).AND.
     &     (CPFEM_mat(CPFEM_mn,2).NE.2)) then
     call CPFEM_error(520)
  endif
    enddo
    CPFEM_tauc_slip_old(:,:,CPFEM_en,CPFEM_in)=s0_slip
 endif

!    ************************************
!    ***   CP-FEM Calculation   *** 
!    ************************************
!    *** Reinitialization of stress and consistent tangent ***
 CPFEM_s=0
 CPFEM_d=0

!    *** Loop over all the components ***
 do iori=1,nori

!    *** Initialization of the matrices for t=t0 ***
 Fp_old=CPFEM_Fp_old(:,:,iori,CPFEM_en,CPFEM_in)
 tauc_slip_old=CPFEM_tauc_slip_old(:,iori,CPFEM_en,CPFEM_in)
 tauc_slip_new=tauc_slip_old
 g_old=CPFEM_g_old(:,iori,CPFEM_en,CPFEM_in)
 Tstar_v=CPFEM_sigma_old(:,iori,CPFEM_en,CPFEM_in)
 p10=CPFEM_ini_ori(1,iori,CPFEM_en,CPFEM_in)
 P0=CPFEM_ini_ori(2,iori,CPFEM_en,CPFEM_in)
 p20=CPFEM_ini_ori(3,iori,CPFEM_en,CPFEM_in)
 vf=CPFEM_mat(CPFEM_mn,7*iori+2)

!    *** Calculation of the solution at t=t1 ***
 if (modulo(CPFEM_cn,ijaco).EQ.0) then
    call CPFEM_stress(CPFEM_tinc,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new,
     &     g_old,g_new,tauc_slip_old, 
     &     tauc_slip_new,
     &     Tstar_v,cs,cd,p11,P1,p21,dgmaxc,1,iconv,ising,
     &     icut,CPFEM_en,CPFEM_in,CPFEM_inc)
!    ***   Evaluation of ising      ***
!    *** ising=2 => singular matrix in jacobi calculation ***
!    ***      => use old jacobi        ***
    if (ising==2) then
  jpara=0
    endif 
!    *** Calculation of the consistent tangent ***
    CPFEM_d=CPFEM_d+vf*cd
 else
    call CPFEM_stress(CPFEM_tinc,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new, 
     &     g_old,g_new,tauc_slip_old,
     &     tauc_slip_new,
     &     Tstar_v,cs,hh6,p11,P1,p21,dgmaxc,0,iconv,
     &     ising,icut,CPFEM_en,CPFEM_in,CPFEM_inc)
    jpara=0
 endif

!    *** Cases of unsuccessful calculations *** 
!    ***   Evaluation od ising  ***
!    *** ising!=0 => singular matrix ***
 if (ising==1) then
    write(6,*) 'Singular matrix!'
    write(6,*) 'Integration point: ',CPFEM_in
    write(6,*) 'Element: ',CPFEM_en
    call CPFEM_error(700)
    CPFEM_timefactor=1.e5_pRe
    return
 endif

!    ***   Evaluation of icut   ***
!    *** icut!=0 => too many cutbacks ***
 if (icut==1) then
    write(6,*) 'Too many cutbacks'
    write(6,*) 'Integration point: ',CPFEM_in
    write(6,*) 'Element: ',CPFEM_en
    call CPFEM_error(600)
    CPFEM_timefactor=1.e5_pRe
    return
 endif

!    ***  Evaluation of iconv ***
!    *** iconv!=0 => no convergence ***
 if (iconv==1) then
    write(6,*) 'Inner loop did not converged!'
    write(6,*) 'Integration point: ',CPFEM_in
    write(6,*) 'Element:',CPFEM_en
    call CPFEM_error(600)
    CPFEM_timefactor=1.e5_pRe
    return
 else 
    if (iconv==2) then
  write(6,*) 'Outer loop did not converged!'
  write(6,*) 'Integration point: ',CPFEM_in
  write(6,*) 'Element: ',CPFEM_en
  call CPFEM_error(600)
  CPFEM_timefactor=1.e5_pRe
  return
    endif
 endif

!    *** Update the differents matrices for t=t1 ***
 CPFEM_Fp_new(:,:,iori,CPFEM_en,CPFEM_in)=Fp_new
 CPFEM_tauc_slip_new(:,iori,CPFEM_en,CPFEM_in)=tauc_slip_new
 CPFEM_g_new(:,iori,CPFEM_en,CPFEM_in)=g_new
 CPFEM_sigma_new(:,iori,CPFEM_en,CPFEM_in)=Tstar_v

!    *** Calculation of the misorientation ***    
 phi1mis(1)=p10
 PHImis(1)=P0
 phi2mis(1)=p20
 phi1mis(2)=p11
 PHImis(2)=P1
 phi2mis(2)=p21
 call CPFEM_misori(phi1mis,PHImis,phi2mis,orimis)

!    *** Update the results plotted in MENTAT ***
 CPFEM_results(1,iori,CPFEM_en,CPFEM_in)=p11
 CPFEM_results(2,iori,CPFEM_en,CPFEM_in)=P1
 CPFEM_results(3,iori,CPFEM_en,CPFEM_in)=p21
 CPFEM_results(4,iori,CPFEM_en,CPFEM_in)=orimis
 CPFEM_results(5,iori,CPFEM_en,CPFEM_in)=sum(g_new)
 CPFEM_results(7,iori,CPFEM_en,CPFEM_in)=sum(tauc_slip_new)/nslip
 CPFEM_results(21,iori,CPFEM_en,CPFEM_in)=vf

!    *** Evaluation of the maximum shear ***
 dgmax=max(dgmax,dgmaxc)
!    *** Evaluation of the average Cauchy stress ***    
 CPFEM_s=CPFEM_s+vf*cs

 enddo
!    *** End of the loop over the components ***
!    *************************************
!    *** End of the CP-FEM Calculation *** 
!    *************************************

!    *** Approximate relative element thickness ***
 call CPFEM_thick(CPFEM_ffn1,CPFEM_en,CPFEM_in)
!    *** Restoration of the old jacobian if necessary ***
 if (jpara==0) then
    CPFEM_d=CPFEM_jaco_old(:,:,CPFEM_en,CPFEM_in)
 else
!    *** Store the new jacobian ***     
    CPFEM_jaco_old(:,:,CPFEM_en,CPFEM_in)=CPFEM_d
 endif    
!    *** Calculate timefactor ***
 CPFEM_timefactor=dgmax/dgs
         
 return
 end


     
 subroutine CPFEM_stress(
     &CPFEM_tinc,
     &CPFEM_ffn,
     &CPFEM_ffn1,
     &Fp_old,
     &Fp_new,
     &g_old,
     &g_new,
     &tauc_slip_old,
     &tauc_slip_new,
     &Tstar_v,
     &cs,
     &dcs_de,
     &phi1,
     &PHI,
     &phi2,
     &dgmaxc,    
     &isjaco, 
     &iconv,
     &ising,
     &icut,
     &CPFEM_en,
     &CPFEM_in,
     &CPFEM_inc
     &)
c********************************************************************
c This routine calculates the stress for a single component
c and manages the independent time incrmentation
c********************************************************************
 use mpie
 use prec, only: pRe,pIn
 implicit none

!    *** Definition of variables ***
 integer(pIn) isjaco,iconv,ising,icut,CPFEM_en,CPFEM_in,CPFEM_inc
 real(pRe) CPFEM_tinc,CPFEM_ffn(3,3),CPFEM_ffn1(3,3),Fp_old(3,3),
     &      Fp_new(3,3),g_old(nslip),g_new(nslip),
     &      tauc_slip_old(nslip),tauc_slip_new(nslip),
     &      Tstar_v(6),
     &      cs(6),dcs_de(6,6),phi1,PHI,phi2,dgmaxc
 integer(pIn) jcut
 real(pRe) Tstar_v_h(6),tauc_slip_new_h(nslip),
     &      dt_i,delta_Fg(3,3),Fg_i(3,3),
     &      tauc_slip_new_i(nslip),time,mm(6,6)

!    *** Numerical parameters ***
 integer(pIn), parameter :: ncut=7_pIn
 icut=0

!    *** First attempt to calculate Tstar and tauc with initial timestep ***
 Tstar_v_h=Tstar_v
 tauc_slip_new_h=tauc_slip_new
 call CPFEM_stress_int(CPFEM_tinc,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new,
     &      g_old,g_new,tauc_slip_old,
     &      tauc_slip_new,
     &      Tstar_v,cs,dcs_de,phi1,PHI,phi2,dgmaxc,
     &      isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc)
 if ((iconv==0).AND.(ising==0)) then
    return
 endif

!    *** Calculation of stress and resistences with a cut timestep ***
!    ***      when first try did not converge    ***  
 jcut=1_pIn
 dt_i=0.5*CPFEM_tinc
 delta_Fg=0.5*(CPFEM_ffn1-CPFEM_ffn)
 Fg_i=CPFEM_ffn+delta_Fg
 Tstar_v=Tstar_v_h
 tauc_slip_new_i=tauc_slip_new_h
!    *** Start time ***
 time=dt_i
 do while (time.LE.CPFEM_tinc)
    call CPFEM_stress_int(time,CPFEM_ffn,Fg_i,Fp_old,Fp_new,g_old,
     &    g_new,tauc_slip_old,
     &    tauc_slip_new_i,
     &    Tstar_v,cs,mm,phi1,PHI,
     &    phi2,dgmaxc,0_pIn,iconv,ising,CPFEM_en,
     &    CPFEM_in,CPFEM_inc)
 if ((iconv==0).AND.(ising==0)) then
   time=time+dt_i
   Fg_i=Fg_i+delta_Fg
   Tstar_v_h=Tstar_v
   tauc_slip_new_h=tauc_slip_new_i
     else
   jcut=jcut+1
   if (jcut.GT.ncut) then
      icut=1
      return
   endif
   dt_i=0.5*dt_i
   time=time-dt_i
   delta_Fg=0.5*delta_Fg
   Fg_i=Fg_i-delta_Fg
   Tstar_v=Tstar_v_h
   tauc_slip_new_i=tauc_slip_new_h
     endif
 enddo

!    *** Final calculation of stress and resistences withb full timestep ***
 tauc_slip_new=tauc_slip_new_i
 call CPFEM_stress_int(CPFEM_tinc,CPFEM_ffn,CPFEM_ffn1,Fp_old,Fp_new,
     &      g_old,g_new,tauc_slip_old,
     &      tauc_slip_new,
     &      Tstar_v,cs,dcs_de,phi1,PHI,phi2,dgmaxc,
     &      isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc) 
 return
 end



 subroutine CPFEM_stress_int(
     &dt,      ! Time increment
     &Fg_old,     ! Old global deformation gradient
     &Fg_new,     ! New global deformation gradient
     &Fp_old,     ! Old plastic deformation gradient
     &Fp_new,     ! New plastic deformation gradient
     &g_old,     ! Old cumulative plastic strain of a slip system
     &g_new,     ! New cumulative plastic strain of a slip system
     &tauc_slip_old,  ! Old resistence of a slip system
     &tauc_slip_new,  ! New resistence of a slip system
     &Tstar_v,     ! Second Piola-Kirschoff stress tensor
     &cs,      ! Cauchy stress vector
     &dcs_de,     ! Consistent tangent
     &phi1,      ! Euler angle phi1
     &PHI,      ! Euler angle PHI
     &phi2,      ! Euler angle phi2
     &dgmaxc,
     &isjaco,
     &iconv,
     &ising,
     &CPFEM_en,
     &CPFEM_in,
     &CPFEM_inc
     &)     
c********************************************************************
c This routine calculates the stress for a single component
c it is based on the paper by Kalidindi et al.:
c J. Mech. Phys, Solids Vol. 40, No. 3, pp. 537-569, 1992
c it is modified to use anisotropic elasticity matrix
c********************************************************************
 use mpie
 use prec
 implicit none

!    *** Definition of variables *** 
 integer(pIn) isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc
 real(pRe) dt,Fg_old(3,3),Fg_new(3,3),Fp_old(3,3),Fp_new(3,3),
     &      g_old(nslip),g_new(nslip),
     &      tauc_slip_old(nslip),tauc_slip_new(nslip),
     &      Tstar_v(6),
     &      cs(6),dcs_de(6,6),phi1,PHI,phi2,dgmaxc
 integer(pIn) ic     
 real(pRe) gdot_slip(nslip),Fe(3,3),R(3,3),
     &      U(3,3),de(3,3),tauc2(nslip),Fp2(3,3),
     &      sgm2(6),cs1(6),dF(3,3),Fg2(3,3),dev(6)    
!    *** Numerical parameters ***
 real(pRe), parameter :: pert_ct=1.0e-5_pRe  
 
!    *** Error treatment ***
 dgmaxc=0
 iconv=0
 ising=0   

!    *********************************************
!    ***  Calculation of the new Cauchy stress ***
!    *********************************************

!    *** Call Newton-Raphson method ***
 call NEWTON_RAPHSON(dt,Fg_old,Fg_new,Fp_old,Fp_new,Fe,gdot_slip,
     &     tauc_slip_old,tauc_slip_new,  
     &     Tstar_v,cs,iconv,ising)
 
!    *** Calculation of the new orientation ***
 call math_pDecomposition(Fe,U,R,ising)
 if (ising==1) then
    return
 endif
 call math_RtoEuler(transpose(R),phi1,PHI,phi2)

!    *** Evaluation of the maximum slip shear ***
 dgmaxc=maxval(abs(gdot_slip*dt))
 g_new=g_old+abs(gdot_slip)*dt

!    *** Choice of the calculation of the consistent tangent ***     
 if (isjaco==0) then
    return
 endif

!    *********************************************
!    *** Calculation of the consistent tangent ***
!    *********************************************

!    *** Calculation of the consistent tangent with perturbation ***
!    *** Perturbation on the component of Fg ***    
 do ic=1,6    

!    *** Method of small perturbation
 dev=0
 if(ic.le.3) dev(ic)=pert_ct
 if(ic.gt.3) dev(ic)=pert_ct/2
 call CPFEM_conv6to33(dev,de)
 dF=matmul(de,Fg_old) 
 Fg2=Fg_new+dF
 sgm2=Tstar_v
 tauc2=tauc_slip_new

!    *** Calculation of the perturbated Cauchy stress ***
 call NEWTON_RAPHSON(dt,Fg_old,Fg2,Fp_old,Fp2,Fe,gdot_slip,
     &     tauc_slip_old,tauc2, 
     &     sgm2,cs1,iconv,ising)     
  
!    *** Consistent tangent ***
 dcs_de(:,ic)=(cs1-cs)/pert_ct     
 enddo
    
 return
 end


 

 subroutine NEWTON_RAPHSON(
     &dt, 
     &Fg_old,     
     &Fg_new,     
     &Fp_old,     
     &Fp_new,
     &Fe,       
     &gdot_slip,    
     &tauc_slip_old,  
     &tauc_slip_new,  
     &Tstar_v, 
     &cs,
     &iconv,
     &ising          
     &) 
!***********************************************************************
!***        NEWTON-RAPHSON Calculation      ***
!***********************************************************************
 use mpie
 use prec
 implicit none

!    *** Definition of variables *** 
 integer(pIn) isjaco,iconv,ising,CPFEM_en,CPFEM_in,CPFEM_inc
 real(pRe) dt,Fg_old(3,3),Fg_new(3,3),Fp_old(3,3),Fp_new(3,3),
     &      g_old(nslip),g_new(nslip),
     &      tauc_slip_old(nslip),tauc_slip_new(nslip),
     &      Tstar_v(6),cs(6),dcs_de(6,6),phi1,PHI,phi2,dgmaxc
 integer(pIn) i,j,k,iouter,iinner,ijac,ic      
 real(pRe) invFp_old(3,3),det,A(3,3),Estar0_v(6),Tstar0_v(6),
     &      mm(3,3),mm1(3,3),vv(6),Dslip(6,nslip),
     &      tau_slip(nslip),gdot_slip(nslip),
     &      R1(6),norm1,Tstar_v_per(6),R1_per(6),
     &      Jacobi(6,6),invJacobi(6,6),dTstar_v(6),R2(nslip),
     &      dtauc_slip(nslip),norm2,dLp(3,3),
     &      Estar(3,3),Estar_v(6),invFp_new(3,3),
     &      invFp2(3,3),Lp(3,3),Fe(3,3),
     &      R(3,3),U(3,3),dgdot_dtaucslip(nslip)
 real(pRe) de(3,3),dev(6),tauc2(nslip),fp2(3,3),
     &      sgm2(6),cs1(6),df(3,3),
     &      fg2(3,3),tauc_old(nslip),crite,tol_in,tol_out

!    *** Numerical parameters ***
 integer(pIn), parameter :: nouter=50
 real(pRe), parameter :: tol_outer=1.0e-4_pRe
 integer(pIn), parameter :: ninner=2000
 real(pRe), parameter :: tol_inner=1.0e-3_pRe   
 real(pRe), parameter :: eta=13.7_pRe 
 integer(pIn), parameter :: numerical=0
 real(pRe), parameter :: pert_nr=1.0e-8_pRe  
 crite=eta*s0_slip/n_slip

!    *** Tolerences ***
 tol_in=tol_inner*s0_slip
 tol_out=tol_outer*s0_slip 

!    *** Error treatment ***
 dgmaxc=0
 iconv=0
 ising=0

!    *** Calculation of Fp_old(-1) ***
 invFp_old=Fp_old   
 call invert(invFp_old,3,0,0,det,3)
 if (det==0.0_pRe) then
    ising=1
    return
 endif

!    *** Calculation of A and T*0 (see Kalidindi) ***
 A=matmul(transpose(matmul(Fg_new,invFp_old)),
     &    matmul(Fg_new,invFp_old))
 call CPFEM_conv33to6((A-I3)/2,Estar0_v)
 Tstar0_v=matmul(Cslip_66,Estar0_v)

!    *** Calculation of Dslip (see Kalidindi) *** 
 do i=1,nslip
    mm=matmul(A,Sslip(i,:,:))
    mm1=(mm+transpose(mm))/2
    vv = math_33to6(mm1)
    Dslip(:,i)=matmul(Cslip_66,vv)
 enddo

!    *** Second level of iterative procedure: Resistences ***
 do iouter=1,nouter
!    *** First level of iterative procedure: Stresses ***
 do iinner=1,ninner

!    *** Calculation of gdot_slip ***    
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
 if (norm1.LT.tol_in) then 
    goto 100
 endif
 
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
    Jacobi(i,i)=1.0_pRe+Jacobi(i,i)
 enddo
 endif
!    *** End of the Jacobi calculation ***

!    *** Inversion of the Jacobi matrix ***
 invJacobi=Jacobi
 call invert(invJacobi,6,0,0,det,6)
 if (det==0.0_pRe) then
    do i=1,6
  Jacobi(i,i)=1.05d0*maxval(Jacobi(i,:))
    enddo
    invJacobi=Jacobi
    call invert(invJacobi,6,0,0,det,6)
    if (det==0.0_pRe) then
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

 enddo
 iconv=1
 return
!    *** End of the first level of iterative procedure ***

 100 continue

 call hardening(tauc_slip_new,gdot_slip,dtauc_slip)

!    *** Arrays of residuals ***
 R2=tauc_slip_new-tauc_slip_old-dtauc_slip*dt 
 norm2=maxval(abs(R2))
 if (norm2.LT.tol_out) then
    goto 200
 endif
 tauc_slip_new=tauc_slip_old+dtauc_slip*dt
 enddo  
 iconv=2
 return
!    *** End of the second level of iterative procedure ***

 200 continue

 call plastic_vel_grad(dt,tau_slip,tauc_slip_new,Lp)
 
!    *** Calculation of Fp(t+dt) (see Kalidindi) ***
 dLp=I3+Lp*dt
 Fp_new=matmul(dLp,Fp_old)
 call CPFEM_determ(Fp_new,det)
 Fp_new=Fp_new/det**(1.0_pRe/3.0_pRe)

!    *** Calculation of F*(t+dt) (see Kalidindi) ***
 invFp_new=Fp_new
 call invert(invFp_new,3,0,0,det,3)
 if (det==0.0_pRe) then
    ising=1
    return
 endif
 Fe=matmul(Fg_new,invFp_new)

!    *** Calculation of Estar ***
 Estar=0.5_pRe*(matmul(transpose(Fe),Fe)-I3)
 call CPFEM_conv33to6(Estar,Estar_v)

!    *** Calculation of the Cauchy stress ***
 call cauchy_stress(Estar_v,Fe,cs)

 return
 end


 
 end module