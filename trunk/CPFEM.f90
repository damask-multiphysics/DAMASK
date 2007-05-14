
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

!*********************************************************
!***    allocate the arrays defined in module CPFEM    ***
!***    and initialize them                            ***
!*********************************************************
 SUBROUTINE CPFEM_init()
!
 use prec, only: pReal,pInt
 use math, only: math_EulertoR
 use mesh
 use constitutive
!
 implicit none

 integer(pInt) e,i,g
!
!    *** mpie.marc parameters ***
 allocate(CPFEM_ffn_all   (3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn_all    = 0.0_pReal
 allocate(CPFEM_ffn1_all  (3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1_all   = 0.0_pReal
 allocate(CPFEM_stress_all(  6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_stress_all = 0.0_pReal
 allocate(CPFEM_jacobi_all(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_jacobi_all = 0.0_pReal
!
!    *** User defined results !!! MISSING incorporate consti_Nresults ***
 allocate(CPFEM_results(CPFEM_Nresults+constitutive_maxNresults,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 CPFEM_results = 0.0_pReal
!
!    *** Second Piola-Kirchoff stress tensor at (t=t0) and (t=t1) ***
 allocate(CPFEM_sigma_old(6,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_sigma_old = 0.0_pReal
 allocate(CPFEM_sigma_new(6,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_sigma_new = 0.0_pReal
!
!    *** Plastic deformation gradient at (t=t0) and (t=t1) ***
 allocate(CPFEM_Fp_old(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems))
 forall (e=1:mesh_NcpElems,i=1:mesh_maxNips,g=1:constitutive_maxNgrains) &
   CPFEM_Fp_old(:,:,g,i,e) = math_EulerToR(constitutive_EulerAngles(:,g,i,e))  ! plastic def gradient reflects init orientation
 allocate(CPFEM_Fp_new(3,3,constitutive_maxNgrains,mesh_maxNips,mesh_NcpElems)) ; CPFEM_Fp_new = 0.0_pReal
!
!    *** Old jacobian (consistent tangent) ***
 allocate(CPFEM_jaco_old(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_jaco_old = 0.0_pReal
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

 END SUBROUTINE
!
!
!***********************************************************************
!***    perform initialization at first call, update variables and   ***
!***    call the actual material model                               ***
!***********************************************************************
 SUBROUTINE CPFEM_general(ffn, ffn1, CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_dt, CPFEM_en, CPFEM_in)
!
 use prec, only: pReal,pInt
 use math, only: math_init
 use mesh, only: mesh_init,mesh_FEasCP
 use constitutive, only: constitutive_init,constitutive_state_old,constitutive_state_new
 implicit none
!
 real(pReal)   ffn(3,3), ffn1(3,3), CPFEM_dt
 integer(pInt) CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_en, CPFEM_in, cp_en
!
! initialization step
 if (CPFEM_first_call) then
! three dimensional stress state ?
    call math_init()
    call mesh_init()
    call constitutive_init()
    call CPFEM_init()
    CPFEM_first_call = .false.
 endif
 if (CPFEM_inc==CPFEM_inc_old) then ! not a new increment
! case of a new subincrement:update starting with subinc 2
     if (CPFEM_subinc > CPFEM_subinc_old) then
        CPFEM_sigma_old        = CPFEM_sigma_new
        CPFEM_Fp_old           = CPFEM_Fp_new
        constitutive_state_old = constitutive_state_new
        CPFEM_subinc_old       = CPFEM_subinc
    endif
 else                               ! new increment
    CPFEM_sigma_old         = CPFEM_sigma_new
    CPFEM_Fp_old            = CPFEM_Fp_new
    constitutive_state_old  = constitutive_state_new
    CPFEM_inc_old           = CPFEM_inc
    CPFEM_subinc_old        = 1_pInt
 endif

 cp_en = mesh_FEasCP('elem',CPFEM_en)
 CPFEM_ffn_all(:,:,CPFEM_in, cp_en)  = ffn
 CPFEM_ffn1_all(:,:,CPFEM_in, cp_en) = ffn1
 call CPFEM_stressIP(CPFEM_cn, CPFEM_dt, cp_en, CPFEM_in)
 return

 END SUBROUTINE


!**********************************************************
!***  calculate the material behaviour at IP level      ***
!**********************************************************
 SUBROUTINE CPFEM_stressIP(&
     CPFEM_cn,&       ! Cycle number
     CPFEM_dt,&       ! Time increment (dt)
     cp_en,&          ! Element number
     CPFEM_in)        ! Integration point number

 use prec, only: pReal,pInt,ijaco,nCutback
 use math, only: math_pDecomposition,math_RtoEuler
 use IO,   only: IO_error
 use mesh, only: mesh_element
 use constitutive
!
 implicit none

 integer(pInt), parameter :: i_now = 1_pInt,i_then = 2_pInt
 character(len=128) msg
 integer(pInt) CPFEM_cn,cp_en,CPFEM_in,grain,i
 logical updateJaco,error
 real(pReal) CPFEM_dt,dt,t,volfrac
 real(pReal), dimension(6) :: cs,Tstar_v
 real(pReal), dimension(6,6) :: cd
 real(pReal), dimension(3,3) :: Fe,U,R,deltaFg
 real(pReal), dimension(3,3,2) :: Fg,Fp
 real(pReal), dimension(constitutive_maxNstatevars,2) :: state

 updateJaco = (mod(CPFEM_cn,ijaco)==0)   ! update consistent tangent every ijaco'th iteration

 CPFEM_stress_all(:,CPFEM_in,cp_en) = 0.0_pReal                  ! average Cauchy stress
 if (updateJaco) CPFEM_jaco_old(:,:,CPFEM_in,cp_en) = 0.0_pReal  ! average consistent tangent

! -------------- grain loop -----------------
 do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
! -------------------------------------------

   i = 0_pInt                         ! cutback counter
   state(:,i_now) = constitutive_state_old(:,grain,CPFEM_in,cp_en)
   Fg(:,:,i_now)  = CPFEM_ffn_all(:,:,CPFEM_in,cp_en)
   Fp(:,:,i_now)  = CPFEM_Fp_old(:,:,grain,CPFEM_in,cp_en)

   deltaFg = CPFEM_ffn1_all(:,:,CPFEM_in,cp_en)-CPFEM_ffn_all(:,:,CPFEM_in,cp_en)
   dt = CPFEM_dt

   Tstar_v = 0.0_pReal                ! fully elastic initial guess
   Fg(:,:,i_then) = Fg(:,:,i_now)
   state(:,i_then) = 0.0_pReal        ! state_old as initial guess
   t = 0.0_pReal

! ------- crystallite integration -----------
   do
! -------------------------------------------
     if (t+dt < CPFEM_dt) then          ! intermediate solution
       t = t+dt                         ! next time inc
       Fg(:,:,i_then) = Fg(:,:,i_then)+deltaFg  ! corresponding Fg
     else                               ! full step solution
       t = CPFEM_dt                     ! final time
       Fg(:,:,i_then) = CPFEM_ffn1_all(:,:,CPFEM_in,cp_en) ! final Fg
     endif

     call CPFEM_stressCrystallite(msg,cs,cd,Tstar_v,Fp(:,:,i_then),Fe,state(:,i_then),&
                                  dt,cp_en,CPFEM_in,grain,updateJaco .and. t==CPFEM_dt,&
                                  Fg(:,:,i_now),Fg(:,:,i_then),Fp(:,:,i_now),state(:,i_now))
     if (msg == 'ok') then             ! solution converged
       if (t == CPFEM_dt) exit         ! reached final "then"
     else                              ! solution not found
       i = i+1_pInt                    ! inc cutback counter
       if (i > nCutback) then          ! limit exceeded?
         write(6,*) 'cutback limit --> '//msg
         write(6,*) 'Grain:             ',grain
         write(6,*) 'Integration point: ',CPFEM_in
         write(6,*) 'Element:           ',mesh_element(1,cp_en)
         call IO_error(600)
         return                        ! byebye
       else
         t = t-dt                        ! rewind time
         Fg(:,:,i_then) = Fg(:,:,i_then)-deltaFg ! rewind Fg
         dt = 0.5_pReal*dt               ! cut time-step in half
         deltaFg = 0.5_pReal*deltaFg     ! cut Fg-step in half
       endif
     endif
   enddo    ! crystallite integration (cutback loop)

! ---- update crystallite matrices at t = t1 ----
   CPFEM_Fp_new(:,:,grain,CPFEM_in,cp_en)         = Fp(:,:,i_then)
   constitutive_state_new(:,grain,CPFEM_in,cp_en) = state(:,i_then)
   CPFEM_sigma_new(:,grain,CPFEM_in,cp_en)        = Tstar_v
! ---- update results plotted in MENTAT ----
   call math_pDecomposition(Fe,U,R,error) ! polar decomposition
   if (error) then
     write(6,*) 'polar decomposition'
     write(6,*) 'Grain:             ',grain
     write(6,*) 'Integration point: ',CPFEM_in
     write(6,*) 'Element:           ',mesh_element(1,cp_en)
     call IO_error(600)
     return
   endif
   CPFEM_results(1:3,grain,CPFEM_in,cp_en) = math_RtoEuler(transpose(R))        ! orientation
   CPFEM_results(4:3+constitutive_Nresults(grain,CPFEM_in,cp_en),grain,CPFEM_in,cp_en) = &
     constitutive_post_results(Tstar_v,state(:,i_then),CPFEM_dt,grain,CPFEM_in,cp_en)

! ---- contribute to IP result ----
    volfrac = constitutive_matVolFrac(grain,CPFEM_in,cp_en)*constitutive_texVolFrac(grain,CPFEM_in,cp_en)
    CPFEM_stress_all(:,CPFEM_in,cp_en) = CPFEM_stress_all(:,CPFEM_in,cp_en)+volfrac*cs                  ! average Cauchy stress
    if (updateJaco) CPFEM_jaco_old(:,:,CPFEM_in,cp_en) = CPFEM_jaco_old(:,:,CPFEM_in,cp_en)+volfrac*cd  ! average consistent tangent

 enddo    ! grain loop

 return
 END SUBROUTINE


!********************************************************************
! Calculates the stress for a single component
! it is based on the paper by Kalidindi et al.:
! J. Mech. Phys, Solids Vol. 40, No. 3, pp. 537-569, 1992
! it is modified to use anisotropic elasticity matrix
!********************************************************************
 subroutine CPFEM_stressCrystallite(&
     msg,&        ! return message
     cs,&         ! Cauchy stress vector
     dcs_de,&     ! consistent tangent
     Tstar_v,&    ! second Piola-Kirchoff stress tensor
     Fp_new,&     ! new plastic deformation gradient
     Fe_new,&     ! new "elastic" deformation gradient
     state_new,&  ! new state variable array
!
     dt,&         ! time increment
     cp_en,&      ! element number
     CPFEM_in,&   ! integration point number
     grain,&      ! grain number
     updateJaco,& ! boolean to calculate Jacobi matrix
     Fg_old,&     ! old global deformation gradient
     Fg_new,&     ! new global deformation gradient
     Fp_old,&     ! old plastic deformation gradient
     state_old)   ! old state variable array

 use prec, only: pReal,pInt,pert_e
 use constitutive, only: constitutive_Nstatevars
 use math, only: math_Mandel6to33,mapMandel
 implicit none

 character(len=*) msg
 logical updateJaco
 integer(pInt) cp_en,CPFEM_in,grain,i
 real(pReal) dt
 real(pReal), dimension(3,3) :: Fg_old,Fg_new,Fg_pert,Fp_old,Fp_new,Fp_pert,Fe_new,Fe_pert,E_pert
 real(pReal), dimension(6)   :: cs,Tstar_v,Tstar_v_pert
 real(pReal), dimension(6,6) :: dcs_de
 real(pReal), dimension(constitutive_Nstatevars(grain,CPFEM_in,cp_en)) :: state_old,state_new,state_pert

 call CPFEM_timeIntegration(msg,Fp_new,Fe_new,Tstar_v,state_new, &   ! def gradients and PK2 at end of time step
                            dt,cp_en,CPFEM_in,grain,Fg_new,Fp_old,state_old)
 if (msg /= 'ok') return
 cs = CPFEM_CauchyStress(Tstar_v,Fe_new)    ! Cauchy stress

 if (updateJaco) then                       ! consistent tangent using numerical perturbation of Fg
   do i = 1,6                               ! Fg component
     E_pert = 0.0_pReal
     E_pert(mapMandel(1,i),mapMandel(2,i)) = E_pert(mapMandel(1,i),mapMandel(2,i)) + pert_e/2.0_pReal
     E_pert(mapMandel(2,i),mapMandel(1,i)) = E_pert(mapMandel(2,i),mapMandel(1,i)) + pert_e/2.0_pReal

     Fg_pert = Fg_new+matmul(E_pert,Fg_old) ! perturbated Fg
     Tstar_v_pert = Tstar_v                 ! initial guess from end of time step
     state_pert = state_new                 ! initial guess from end of time step

     call CPFEM_timeIntegration(msg,Fp_pert,Fe_pert,Tstar_v_pert,state_pert, &
                                dt,cp_en,CPFEM_in,grain,Fg_pert,Fp_old,state_old)
     if (msg /= 'ok') then
       msg = 'consistent tangent --> '//msg
       return
     endif
! Remark: (perturbated) Cauchy stress is Mandel hence dcs_de(:,4:6) is too large by sqrt(2)
     dcs_de(:,i) = (CPFEM_CauchyStress(Tstar_v_pert,Fe_pert)-cs)/pert_e
   enddo
 endif

 return

 END SUBROUTINE


!***********************************************************************
!***     fully-implicit two-level time integration                   ***
!***********************************************************************
 SUBROUTINE CPFEM_timeIntegration(&
     msg,&              ! return message
     Fp_new,&           ! new plastic deformation gradient
     Fe_new,&           ! new "elastic" deformation gradient
     Tstar_v,&          ! 2nd PK stress (taken as initial guess if /= 0)
     state_new,&        ! current microstructure at end of time inc (taken as guess if /= 0)
!
     dt,&               ! time increment
     cp_en,&            ! element number
     CPFEM_in,&         ! integration point number
     grain,&            ! grain number
     Fg_new,&           ! new total def gradient
     Fp_old,&           ! former plastic def gradient
     state_old)         ! former microstructure

 use prec, only: pReal,pInt, nState,tol_State,nStress,tol_Stress, crite, nReg
 use constitutive, only: constitutive_Nstatevars,&
                         constitutive_homogenizedC,constitutive_dotState,constitutive_LpAndItsTangent
 use math
 implicit none

 character(len=*) msg
 integer(pInt) cp_en, CPFEM_in, grain
 integer(pInt) iState,iStress,dummy, i,j,k,l,m
 real(pReal) dt,det
 real(pReal), dimension(6) :: Tstar_v,dTstar_v,Rstress
 real(pReal), dimension(6,6) :: C_66,Jacobi,invJacobi
 real(pReal), dimension(3,3) :: Fg_new,Fp_old,Fp_new,Fe_new,invFp_old,invFp_new,Lp,A,B,AB
 real(pReal), dimension(3,3,3,3) :: dLp, LTL
 real(pReal), dimension(constitutive_Nstatevars(grain, CPFEM_in, cp_en)) :: state_old,state_new,dstate,Rstate,RstateS
 logical failed

 msg = 'ok'  ! error-free so far

 call math_invert3x3(Fp_old,invFp_old,det,failed) ! inversion of Fp
 if (failed) then
    msg = 'inversion Fp_old'
    return
 endif

 C_66 = constitutive_HomogenizedC(grain, CPFEM_in, cp_en)
 A = matmul(Fg_new,invFp_old)  ! actually Fe
 A = matmul(transpose(A), A)

! former state guessed, if none specified
 if (all(state_new == 0.0_pReal)) state_new = state_old
 RstateS = state_new
 iState = 0_pInt
! fully elastic guess (Lp = 0), if none specified
 if (all(Tstar_v == 0.0_pReal)) Tstar_v = 0.5_pReal*matmul(C_66,math_Mandel33to6(A-math_I3))
! QUESTION follow former plastic slope to guess better?
 Rstress = Tstar_v

state: do                ! outer iteration: state
         iState = iState+1
         if (iState > nState) then
           msg = 'limit state iteration'
           return
         endif
         iStress = 0_pInt
stress:  do              ! inner iteration: stress
           iStress = iStress+1
           if (iStress > nStress) then      ! too many loops required
             msg = 'limit stress iteration'
             return
           endif
           call constitutive_LpAndItsTangent(Lp,dLp, Tstar_v,state_new,grain,CPFEM_in,cp_en)
           B = math_I3-dt*Lp
           AB = matmul(A,B)
           Rstress = Tstar_v - 0.5_pReal*matmul(C_66,math_Mandel33to6(matmul(transpose(B),AB)-math_I3))
           if (maxval(abs(Tstar_v)) == 0.0_pReal .or. maxval(abs(Rstress/maxval(abs(Tstar_v)))) < tol_Stress) exit stress

!   update stress guess using inverse of dRes/dTstar (Newton--Raphson)
           LTL = 0.0_pReal
           do i=1,3
             do j=1,3
               do k=1,3
                 do l=1,3
                   do m=1,3
!                    LTL(i,j,k,l) = LTL(i,j,k,l) + AB(i,m)*dLp(m,j,k,l) + AB(j,m)*dLp(m,i,l,k)  ! old
                     LTL(i,j,k,l) = LTL(i,j,k,l) + dLp(j,i,m,k)*AB(m,l) + AB(m,i)*dLp(m,j,k,l)   ! new (and correct??)
                   enddo
                 enddo
               enddo
             enddo
           enddo

           Jacobi = math_identity2nd(6) + 0.5_pReal*dt*matmul(C_66,math_Mandel3333to66(LTL))
           j = 0_pInt ; failed = .true.
           do while (failed .and. j <= nReg)
             call math_invert6x6(Jacobi,invJacobi,dummy,failed)
             forall (i=1:6) Jacobi(i,i) = 1.05_pReal*maxval(Jacobi(i,:)) ! regularization
             j = j+1
           enddo
           if (failed) then
             msg = 'regularization Jacobi'
             return
           endif

           dTstar_v = matmul(invJacobi,Rstress)  ! correction to Tstar
           forall(i=1:6, abs(dTstar_v(i)) > crite*maxval(abs(Tstar_v))) &
             dTstar_v(i) = sign(crite*maxval(abs(Tstar_v)),dTstar_v(i))   ! cap to maximum correction
           Tstar_v = Tstar_v-dTstar_v

    enddo stress

    dstate = dt*constitutive_dotState(Tstar_v,state_new,grain,CPFEM_in,cp_en) ! evolution of microstructure
    Rstate = state_new - (state_old+dstate)
    RstateS = 0.0_pReal
    forall (i=1:constitutive_Nstatevars(grain,CPFEM_in,cp_en), state_new(i)/=0.0_pReal) &
      RstateS(i) = Rstate(i)/state_new(i)
    if (maxval(abs(RstateS)) < tol_State) exit state
    state_new = state_old+dstate

 enddo state

 invFp_new = matmul(invFp_old,B)
 call math_invert3x3(invFp_new,Fp_new,det,failed)
 if (failed) then
    msg = 'inversion Fp_new'
    return
 endif
 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal) ! det = det(InvFp_new) !!
 Fe_new = matmul(Fg_new,invFp_new)

 return
 END SUBROUTINE


 FUNCTION CPFEM_CauchyStress(PK_v,Fe)
!***********************************************************************
!***        Cauchy stress calculation                               ***
!***********************************************************************
 use prec, only: pReal,pInt
 use math, only: math_Mandel33to6,math_Mandel6to33,math_det3x3
 implicit none
!    *** Subroutine parameters ***
 real(pReal) PK_v(6), Fe(3,3), CPFEM_CauchyStress(6)

 CPFEM_CauchyStress = math_Mandel33to6(matmul(matmul(Fe,math_Mandel6to33(PK_v)),transpose(Fe))/math_det3x3(Fe))
 return
 END FUNCTION


 END MODULE
