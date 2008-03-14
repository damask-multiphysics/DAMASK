
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
 real(pReal), dimension (:,:),       allocatable :: CPFEM_Temperature
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
 real(pReal), dimension (:,:,:,:),   allocatable :: CPFEM_jacobian
 real(pReal), parameter :: CPFEM_odd_stress = 1e15_pReal, CPFEM_odd_jacobian = 1e50_pReal
 integer(pInt) :: CPFEM_inc_old    = 0_pInt
 integer(pInt) :: CPFEM_subinc_old = 1_pInt
 integer(pInt) :: CPFEM_cycle_old = -1_pInt
 integer(pInt) :: CPFEM_Nresults   = 4_pInt    ! three Euler angles plus volume fraction
 logical :: CPFEM_first_call = .true.

 CONTAINS

!*********************************************************
!***    allocate the arrays defined in module CPFEM    ***
!***    and initialize them                            ***
!*********************************************************
 SUBROUTINE CPFEM_init()
!
 use prec
 use math, only: math_EulertoR, math_I3, math_identity2nd
 use mesh
 use constitutive
!
 implicit none

 integer(pInt) e,i,g
!
!    *** mpie.marc parameters ***
 allocate(CPFEM_Temperature   (mesh_maxNips,mesh_NcpElems)) ; CPFEM_Temperature = 0.0_pReal
 allocate(CPFEM_ffn_all   (3,3,mesh_maxNips,mesh_NcpElems))
 forall(e=1:mesh_NcpElems,i=1:mesh_maxNips) CPFEM_ffn_all(:,:,i,e)             = math_I3
 allocate(CPFEM_ffn1_all  (3,3,mesh_maxNips,mesh_NcpElems)) ; CPFEM_ffn1_all   = CPFEM_ffn_all
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
!    *** FEM jacobian (consistent tangent) ***
 allocate(CPFEM_jacobian(6,6,mesh_maxNips,mesh_NcpElems)) ; CPFEM_jacobian = 0.0_pReal
!
!
!    *** Output to MARC output file ***
 write(6,*)
 write(6,*) 'Arrays allocated:'
 write(6,*) 'CPFEM_Temperature:   ', shape(CPFEM_Temperature)
 write(6,*) 'CPFEM_ffn_all:       ', shape(CPFEM_ffn_all)
 write(6,*) 'CPFEM_ffn1_all:      ', shape(CPFEM_ffn1_all)
 write(6,*) 'CPFEM_stress_all:    ', shape(CPFEM_stress_all)
 write(6,*) 'CPFEM_jacobi_all:    ', shape(CPFEM_jacobi_all)
 write(6,*) 'CPFEM_results:       ', shape(CPFEM_results)
 write(6,*) 'CPFEM_sigma_old:     ', shape(CPFEM_sigma_old)
 write(6,*) 'CPFEM_sigma_new:     ', shape(CPFEM_sigma_new)
 write(6,*) 'CPFEM_Fp_old:        ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:        ', shape(CPFEM_Fp_new)
 write(6,*) 'CPFEM_jacobian:      ', shape(CPFEM_jacobian)
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
 SUBROUTINE CPFEM_general(ffn, ffn1, Temperature, CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_stress_recovery, CPFEM_dt,&
                          CPFEM_en, CPFEM_in, CPFEM_stress, CPFEM_jaco, CPFEM_ngens)
!
 use prec, only: pReal,pInt
 use debug
 use math, only: math_init, invnrmMandel, math_identity2nd, math_Mandel3333to66,math_Mandel33to6,math_Mandel6to33
 use mesh, only: mesh_init,mesh_FEasCP, mesh_NcpElems, FE_Nips, FE_mapElemtype, mesh_element
 use crystal, only: crystal_Init
 use constitutive, only: constitutive_init,constitutive_state_old,constitutive_state_new,material_Cslip_66
 implicit none

 integer(pInt) CPFEM_inc, CPFEM_subinc, CPFEM_cn, CPFEM_en, CPFEM_in, cp_en, CPFEM_ngens, i, e
 real(pReal)   ffn(3,3),ffn1(3,3),Temperature,CPFEM_dt,CPFEM_stress(CPFEM_ngens),CPFEM_jaco(CPFEM_ngens,CPFEM_ngens)
 logical CPFEM_stress_recovery
 
! calculate only every second cycle

 if (mod(CPFEM_cn,2) /= 0) then ! odd cycle: record data for use in even cycle and return stiff result for this odd cycle
    cp_en = mesh_FEasCP('elem',CPFEM_en)
    CPFEM_Temperature(CPFEM_in, cp_en)  = Temperature
    CPFEM_ffn_all(:,:,CPFEM_in, cp_en)  = ffn
    CPFEM_ffn1_all(:,:,CPFEM_in, cp_en) = ffn1
    CPFEM_stress(1:CPFEM_ngens) = CPFEM_odd_stress
    CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
    CPFEM_cycle_old = CPFEM_cn

 else  ! even cycle: really calculate only in first call of new cycle and when in stress recovery

    if (CPFEM_cn /= CPFEM_cycle_old .and. CPFEM_stress_recovery) then
        if (CPFEM_first_call) then             ! initialization step
                                               ! three dimensional stress state ?
            call math_init()
            call mesh_init()
            call crystal_Init()
            call constitutive_init()
            call CPFEM_init()			
            CPFEM_Temperature  = Temperature
            CPFEM_first_call = .false.
            
        endif

        if (CPFEM_inc == CPFEM_inc_old) then   ! not a new increment
            if (CPFEM_subinc > CPFEM_subinc_old) then  ! new subincrement: update starting with subinc 2
                CPFEM_sigma_old        = CPFEM_sigma_new
                CPFEM_Fp_old           = CPFEM_Fp_new
                constitutive_state_old = constitutive_state_new
                CPFEM_subinc_old       = CPFEM_subinc
            endif
        else                                   ! new increment
            CPFEM_sigma_old         = CPFEM_sigma_new
            CPFEM_Fp_old            = CPFEM_Fp_new
            constitutive_state_old  = constitutive_state_new
            CPFEM_inc_old           = CPFEM_inc
            CPFEM_subinc_old        = 1_pInt
        endif
        CPFEM_cycle_old = CPFEM_cn

        debug_cutbackDistribution = 0_pInt     ! initialize debugging data
        debug_InnerLoopDistribution = 0_pInt
        debug_OuterLoopDistribution = 0_pInt

! this shall be done in a parallel loop in the future

        do e=1,mesh_NcpElems
            do i=1,FE_Nips(FE_mapElemtype(mesh_element(2,e)))
                debugger = (e==1 .and. i==1)
                call CPFEM_stressIP(CPFEM_cn, CPFEM_dt, i, e)
            enddo
        enddo

        call debug_info()        ! output of debugging/performance statistics
    end if

! return stress and jacobi
    cp_en = mesh_FEasCP('elem', CPFEM_en)
    CPFEM_stress(1:CPFEM_ngens) = CPFEM_stress_all(1:CPFEM_ngens, CPFEM_in, cp_en)
    CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_jacobian(1:CPFEM_ngens,1:CPFEM_ngens, CPFEM_in, cp_en)
 end if
 return

 END SUBROUTINE


!**********************************************************
!***  calculate the material behaviour at IP level      ***
!**********************************************************
 SUBROUTINE CPFEM_stressIP(&
     CPFEM_cn,&       ! Cycle number
     CPFEM_dt,&       ! Time increment (dt)
     CPFEM_in,&       ! Integration point number
     cp_en)           ! Element number

 use prec, only: pReal,pInt,ijaco,nCutback
 use debug
 use math, only: math_pDecomposition,math_RtoEuler, inDeg, math_I3, math_invert3x3
 use IO,   only: IO_error
 use mesh, only: mesh_element
 use constitutive

 implicit none

 integer(pInt), parameter :: i_now = 1_pInt,i_then = 2_pInt
 character(len=128) msg
 integer(pInt) CPFEM_cn,cp_en,CPFEM_in,grain,i,max_cutbacks
 logical updateJaco,error,cutback
 real(pReal) CPFEM_dt,dt,t,volfrac,det
 real(pReal), dimension(6) :: cs,Tstar_v
 real(pReal), dimension(6,6) :: cd
 real(pReal), dimension(3,3) :: Fe,U,R,deltaFg,invFgthen,invFpnow,Lp
 real(pReal), dimension(3,3,2) :: Fg,Fp
 real(pReal), dimension(constitutive_maxNstatevars,2) :: state

 updateJaco = (mod(CPFEM_cn,2_pInt*ijaco)==0)   ! update consistent tangent every ijaco'th iteration

 CPFEM_stress_all(:,CPFEM_in,cp_en) = 0.0_pReal                  ! average Cauchy stress
 if (updateJaco) CPFEM_jacobian(:,:,CPFEM_in,cp_en) = 0.0_pReal  ! average consistent tangent

! -------------- grain loop -----------------
 do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
! -------------------------------------------

   i = 0_pInt                         ! cutback counter
   max_cutbacks = 0_pInt              ! maximum depth of cut backing
   dt = CPFEM_dt
   state(:,i_now) = constitutive_state_old(:,grain,CPFEM_in,cp_en)
   Fg(:,:,i_now)  = CPFEM_ffn_all(:,:,CPFEM_in,cp_en)
   Fp(:,:,i_now)  = CPFEM_Fp_old(:,:,grain,CPFEM_in,cp_en)
   invFgthen = 0.0_pReal
   invFpnow = 0.0_pReal
   call math_invert3x3(CPFEM_ffn1_all(:,:,CPFEM_in,cp_en),invFgthen,det,error)
   call math_invert3x3(Fp(:,:,i_now),invFpnow,det,error)
   if (dt /= 0.0_pReal) then
     Lp = (math_I3-matmul(Fp(:,:,i_now),matmul(invFgthen,matmul(Fg(:,:,i_now),invFpnow))))/dt ! fully plastic initial guess
   else
     Lp = 0.0_pReal                                          ! fully elastic guess 
   endif

   deltaFg = CPFEM_ffn1_all(:,:,CPFEM_in,cp_en)-CPFEM_ffn_all(:,:,CPFEM_in,cp_en)

   Tstar_v = CPFEM_sigma_old(:,grain,CPFEM_in,cp_en)         ! use last result as initial guess
   Fg(:,:,i_then) = Fg(:,:,i_now)
   Fp(:,:,i_then) = Fp(:,:,i_now)
   state(:,i_then) = 0.0_pReal                               ! state_old as initial guess
   t = 0.0_pReal
   cutback = .false.                                         ! no cutback has happened so far
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

     call CPFEM_stressCrystallite(msg,cs,cd,Tstar_v,Lp,Fp(:,:,i_then),Fe,state(:,i_then),&
                                  t,cp_en,CPFEM_in,grain,updateJaco .and. t==CPFEM_dt,&
                                  Fg(:,:,i_then),Fp(:,:,i_now),state(:,i_now))
     if (msg == 'ok') then             ! solution converged
       if (t == CPFEM_dt) then
	     debug_cutbackDistribution(max_cutbacks+1) = debug_cutbackDistribution(max_cutbacks+1)+1
	     exit                          ! reached final "then"
	   endif
       if (cutback == .false.) then    ! stable solution at current speed?
         dt = 2.0_pReal*dt             ! double time-step
         i = i-1_pInt                  ! dec cutback counter
       endif
       cutback = .false.               ! solution in next step does not derive from a cutback
     else                              ! solution not found
       i = i+1_pInt                    ! inc cutback counter
       max_cutbacks = max(i,max_cutbacks)
       cutback = .true.
       if (i > nCutback) then          ! limit exceeded?
	     debug_cutbackDistribution(nCutback+1) = debug_cutbackDistribution(nCutback+1)+1
         write(6,'(x,a,x,i6,x,a,x,i2,x,a,x,i2)') 'element:',cp_en,'IP:',CPFEM_in,'grain:',grain
         write(6,*) 'cutback limit --> '//msg
         call IO_error(600)
         return                          ! byebye
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
! ---- contribute to IP result ----
    volfrac = constitutive_matVolFrac(grain,CPFEM_in,cp_en)*constitutive_texVolFrac(grain,CPFEM_in,cp_en)
    CPFEM_stress_all(:,CPFEM_in,cp_en) = CPFEM_stress_all(:,CPFEM_in,cp_en)+volfrac*cs                  ! average Cauchy stress
    if (updateJaco) CPFEM_jacobian(:,:,CPFEM_in,cp_en) = CPFEM_jacobian(:,:,CPFEM_in,cp_en)+volfrac*cd  ! average consistent tangent
! ---- update results plotted in MENTAT ----
   call math_pDecomposition(Fe,U,R,error) ! polar decomposition
   if (error) then
     write(6,*) 'polar decomposition'
     write(6,*) 'Grain:             ',grain
     write(6,*) 'Integration point: ',CPFEM_in
     write(6,*) 'Element:           ',mesh_element(1,cp_en)
     call IO_error(650)
     return
   endif
   CPFEM_results(1:3,grain,CPFEM_in,cp_en) = math_RtoEuler(transpose(R))*inDeg        ! orientation
   CPFEM_results(4  ,grain,CPFEM_in,cp_en) = volfrac                                  ! volume fraction of orientation
   CPFEM_results(5:4+constitutive_Nresults(grain,CPFEM_in,cp_en),grain,CPFEM_in,cp_en) = &
     constitutive_post_results(Tstar_v,state(:,i_then),CPFEM_dt,CPFEM_Temperature(CPFEM_in,cp_en),grain,CPFEM_in,cp_en)

 enddo    ! grain loop

 return
 END SUBROUTINE


!********************************************************************
! Calculates the stress for a single component
!********************************************************************
 subroutine CPFEM_stressCrystallite(&
     msg,&        ! return message
     cs,&         ! Cauchy stress vector
     dcs_de,&     ! consistent tangent
     Tstar_v,&    ! second Piola-Kirchhoff stress tensor
     Lp,&         ! guess of plastic velocity gradient
     Fp_new,&     ! new plastic deformation gradient
     Fe_new,&     ! new "elastic" deformation gradient
     state_new,&  ! new state variable array
!
     dt,&         ! time increment
     cp_en,&      ! element number
     CPFEM_in,&   ! integration point number
     grain,&      ! grain number
     updateJaco,& ! boolean to calculate Jacobi matrix
     Fg_new,&     ! new global deformation gradient
     Fp_old,&     ! old plastic deformation gradient
     state_old)   ! old state variable array

 use prec, only: pReal,pInt,pert_Fg
 use debug
 use constitutive, only: constitutive_Nstatevars
 use mesh, only: mesh_element
 use math, only: math_Mandel6to33,math_Mandel33to6,math_Mandel3333to66,&
                 math_I3,math_det3x3,math_invert3x3
 implicit none

 character(len=*) msg
 logical updateJaco,error
 integer(pInt) cp_en,CPFEM_in,grain,i,j,k,l,m,n
 real(pReal) dt,invJ,det
 real(pReal), dimension(3,3,3,3) :: A,H
 real(pReal), dimension(3,3) :: Lp,Lp_pert,Fg_new,Fg_pert,Fp_old,Fp_new,invFp_new,Fp_pert,invFp_pert
 real(pReal), dimension(3,3) :: Fe_new,Fe_pert,Tstar,tau,P,P_pert
 real(pReal), dimension(6)   :: cs,Tstar_v,Tstar_v_pert
 real(pReal), dimension(6,6) :: dcs_de
 real(pReal), dimension(constitutive_Nstatevars(grain,CPFEM_in,cp_en)) :: state_old,state_new,state_pert

 call CPFEM_timeIntegration(msg,Lp,Fp_new,Fe_new,Tstar_v,state_new, &   ! def gradients and PK2 at end of time step
                            dt,cp_en,CPFEM_in,grain,Fg_new,Fp_old,state_old)
 if (msg /= 'ok') return                    ! solution not reached --> report back
 Tstar = math_Mandel6to33(Tstar_v)          ! second PK in intermediate
 tau = matmul(Fe_new,matmul(Tstar,transpose(Fe_new))) ! Kirchhoff stress
 invJ = 1.0_pReal/math_det3x3(Fe_new)       ! inverse dilatation of Fe
 cs = math_Mandel33to6(invJ*tau)            ! Cauchy stress
 if (updateJaco) then                       ! consistent tangent using
                                            ! numerical perturbation of Fg (D. Tjahjanto Diss p.106)
   call math_invert3x3(Fp_new,invFp_new,det,error)
   if (error) then
     msg = 'inversion of Fp_new'
     return
   endif
   P = matmul(Fe_new,&
       matmul(Tstar,transpose(invFp_new)))    ! first PK at center
   do k=1,3
     do l=1,3
       Fg_pert = Fg_new                       ! initialize perturbed Fg
       Fg_pert(k,l) = Fg_pert(k,l) + pert_Fg  ! perturb single component
       Lp_pert    = Lp
       state_pert = state_new                 ! initial guess from end of time step
       call CPFEM_timeIntegration(msg,Lp_pert,Fp_pert,Fe_pert,Tstar_v_pert,state_pert, &
                                  dt,cp_en,CPFEM_in,grain,Fg_pert,Fp_old,state_old)
       if (msg /= 'ok') then
         msg = 'consistent tangent --> '//msg
         return
       endif
     
       call math_invert3x3(Fp_pert,invFp_pert,det,error)
       if (error) then
         msg = 'inversion of Fp_pert'
         return
       endif
       P_pert = matmul(Fe_pert,&
                matmul(math_mandel6to33(Tstar_v_pert),transpose(invFp_pert))) ! perturbed first PK
       A(:,:,k,l) = (P_pert-P)/pert_Fg        ! dP_ij/dFg_kl
     enddo
   enddo
   
   H = 0.0_pReal
   forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
     H(i,j,k,l) = H(i,j,k,l) + &
                  (Fg_new(j,m)*Fg_new(l,n)*A(i,m,k,n) - math_I3(j,l)*Fg_new(i,m)*P(k,m)) + &
                  0.5_pReal*(math_I3(i,k)*tau(j,l) + math_I3(j,l)*tau(i,k) + &
                             math_I3(i,l)*tau(j,k) + math_I3(j,k)*tau(i,l))
   dcs_de = math_Mandel3333to66(invJ*H)     ! Mandel version of stiffness tensor
 endif

 return

 END SUBROUTINE


!***********************************************************************
!***     fully-implicit two-level time integration                   ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 SUBROUTINE CPFEM_timeIntegration(&
     msg,&              ! return message
     Lpguess,&          ! guess of plastic velocity gradient
     Fp_new,&           ! new plastic deformation gradient
     Fe_new,&           ! new "elastic" deformation gradient
     Tstar_v,&          ! 2nd PK stress (taken as initial guess if /= 0)
     state,&            ! current microstructure at end of time inc (taken as guess if /= 0)
!
     dt,&               ! time increment
     cp_en,&            ! element number
     CPFEM_in,&         ! integration point number
     grain,&            ! grain number
     Fg_new,&           ! new total def gradient
     Fp_old,&           ! former plastic def gradient
     state_old)         ! former microstructure

 use prec
 use debug
 use mesh, only: mesh_element
 use constitutive, only: constitutive_Nstatevars,&
                         constitutive_homogenizedC,constitutive_dotState,constitutive_LpAndItsTangent,&
 						 constitutive_Microstructure
 use math
 implicit none

 character(len=*) msg
 integer(pInt) cp_en, CPFEM_in, grain
 integer(pInt) iOuter,iInner,dummy, i,j,k,l,m,n
 real(pReal) dt, det, p_hydro, leapfrog,maxleap
 real(pReal), dimension(6) :: Tstar_v
 
 real(pReal), dimension(9,9) :: dLp,dTdLp,dRdLp,invdRdLp,eye2
 real(pReal), dimension(6,6) :: C_66
 real(pReal), dimension(3,3) :: Fg_new,invFg_new,Fp_new,invFp_new,Fp_old,invFp_old,Fe_new,Fe_old
 real(pReal), dimension(3,3) :: Tstar
 real(pReal), dimension(3,3) :: Lp,Lpguess,Lpguess_old,Rinner,Rinner_old,A,B,BT,AB,BTA
 real(pReal), dimension(3,3,3,3) :: C
 real(pReal), dimension(constitutive_Nstatevars(grain, CPFEM_in, cp_en)) :: state_old,state,ROuter
 logical failed

 msg = 'ok'  ! error-free so far
 
 eye2 = math_identity2nd(9)
 call math_invert3x3(Fp_old,invFp_old,det,failed) ! inversion of Fp_old
 if (failed) then
    msg = 'inversion Fp_old'
    return
 endif
 call math_invert3x3(Fg_new,invFg_new,det,failed) ! inversion of Fg_new
 if (failed) then
    msg = 'inversion Fg_new'
    return
 endif

 Fe_old = matmul(Fg_new,invFp_old)
 A = matmul(transpose(Fe_old), Fe_old)

 if (all(state == 0.0_pReal)) state = state_old    ! former state guessed, if none specified
 iOuter = 0_pInt                                   ! outer counter

Outer: do                ! outer iteration: State
         iOuter = iOuter+1
         if (iOuter > nOuter) then
           msg = 'limit Outer iteration'
		   debug_OuterLoopDistribution(nOuter) = debug_OuterLoopDistribution(nOuter)+1
           return
         endif
		 call constitutive_Microstructure(state,CPFEM_Temperature(CPFEM_in,cp_en),grain,CPFEM_in,cp_en)
		 C_66 = constitutive_HomogenizedC(state, grain, CPFEM_in, cp_en)
         C = math_Mandel66to3333(C_66)       ! 4th rank elasticity tensor
         
         iInner = 0_pInt
         leapfrog = 1.0_pReal                ! correction as suggested by invdRdLp-step
         maxleap = 1024.0_pReal              ! preassign maximum acceleration level

Inner:  do              ! inner iteration: Lp
           iInner = iInner+1
           if (iInner > nInner) then         ! too many loops required
             msg = 'limit Inner iteration'
		     debug_InnerLoopDistribution(nInner) = debug_InnerLoopDistribution(nInner)+1
             return
           endif
           B = math_i3 - dt*Lpguess
           BT = transpose(B)
           AB = matmul(A,B)
           BTA = matmul(BT,A)
           Tstar_v = 0.5_pReal*matmul(C_66,math_mandel33to6(matmul(BT,AB)-math_I3))
           Tstar = math_Mandel6to33(Tstar_v)
           p_hydro=(Tstar_v(1)+Tstar_v(2)+Tstar_v(3))/3.0_pReal
           forall(i=1:3) Tstar_v(i) = Tstar_v(i)-p_hydro       ! subtract hydrostatic pressure
           call constitutive_LpAndItsTangent(Lp,dLp, &
                                             Tstar_v,state,CPFEM_Temperature(CPFEM_in,cp_en),grain,CPFEM_in,cp_en)
           Rinner = Lpguess - Lp                   ! update current residuum
           if (( maxval(abs(Rinner)) < abstol_Inner ) .or. &
               ( any(abs(dt*Lpguess) > relevantStrain) .and. &
                maxval(abs(Rinner/Lpguess),abs(dt*Lpguess) > relevantStrain) < reltol_Inner )&
              ) exit Inner

           ! check for acceleration/deceleration in Newton--Raphson correction
           
           if (leapfrog > 1.0_pReal .and. &
               (sum(Rinner*Rinner) > sum(Rinner_old*Rinner_old) .or. &  ! worse residuum
               sum(Rinner*Rinner_old) < 0.0_pReal)) then                ! residuum changed sign (overshoot)

             maxleap = 0.5_pReal * leapfrog                             ! limit next acceleration
             leapfrog = 1.0_pReal                                       ! grinding halt

           else                                                         ! better residuum

             dTdLp = 0.0_pReal                                          ! calc dT/dLp
             forall (i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
               dTdLp(3*(i-1)+j,3*(k-1)+l) = dTdLp(3*(i-1)+j,3*(k-1)+l) + &
               C(i,j,l,n)*AB(k,n)+C(i,j,m,l)*BTA(m,k)
             dTdLp = -0.5_pReal*dt*dTdLp
           
             dRdLp = eye2 - matmul(dLp,dTdLp)                           ! calc dR/dLp

             invdRdLp = 0.0_pReal
             call math_invert(9,dRdLp,invdRdLp,dummy,failed)            ! invert dR/dLp --> dLp/dR
             if (failed) then
               msg = 'inversion dR/dLp'
               return
             endif

             Rinner_old = Rinner                                        ! remember current residuum
             Lpguess_old = Lpguess                                      ! remember current Lp guess
             if (iInner > 1 .and. leapfrog < maxleap) &
               leapfrog = 2.0_pReal * leapfrog                          ! accelerate
           endif

           Lpguess = Lpguess_old                                        ! start from current guess
           Rinner  = Rinner_old                                         ! use current residuum
           forall (i=1:3,j=1:3,k=1:3,l=1:3) &                           ! leapfrog to updated Lpguess 
             Lpguess(i,j) = Lpguess(i,j) - leapfrog*invdRdLp(3*(i-1)+j,3*(k-1)+l)*Rinner(k,l)
  		   
     enddo Inner
 
       debug_InnerLoopDistribution(iInner) = debug_InnerLoopDistribution(iInner)+1
	   ROuter = state - state_old - &
	            dt*constitutive_dotState(Tstar_v,state,CPFEM_Temperature(CPFEM_in,cp_en),&
	                                     grain,CPFEM_in,cp_en)          ! residuum from evolution of microstructure
       state = state - ROuter                                           ! update of microstructure
       if (maxval(abs(Router/state),state /= 0.0_pReal) < reltol_Outer) exit Outer

   enddo Outer

 debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
 invFp_new = matmul(invFp_old,B)
 call math_invert3x3(invFp_new,Fp_new,det,failed)
 if (failed) then
    msg = 'inversion Fp_new'
    return
 endif
 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal)     ! regularize Fp by det = det(InvFp_new) !!
 Fe_new = matmul(Fg_new,invFp_new)              ! calc resulting Fe
 forall (i=1:3) Tstar_v(i) = Tstar_v(i)+p_hydro ! add hydrostatic component back

 return
 
 END SUBROUTINE


 END MODULE
