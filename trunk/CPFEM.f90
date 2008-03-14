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
 real(pReal), dimension (:,:,:,:),    allocatable :: CPFEM_results
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
!    *** Output to MARC output file ***
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
 write(6,*) 'CPFEM_results:       ', shape(CPFEM_results)
 write(6,*) 'CPFEM_Fp_old:        ', shape(CPFEM_Fp_old)
 write(6,*) 'CPFEM_Fp_new:        ', shape(CPFEM_Fp_new)
 write(6,*)
 call flush(6)
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
 use crystal, only: crystal_Init
 use constitutive, only: constitutive_init,constitutive_state_old,constitutive_state_new,material_Cslip_66
 implicit none
!
 integer(pInt) CPFEM_en, CPFEM_in, cp_en, CPFEM_ngens, i,j,k,l,m,n, e
 real(pReal), dimension (3,3)        :: ffn,ffn1,Kirchhoff_bar
 real(pReal), dimension (3,3,3,3)    :: H_bar
 real(pReal), dimension(CPFEM_ngens) :: CPFEM_stress
 real(pReal), dimension(CPFEM_ngens,CPFEM_ngens) :: CPFEM_jaco
 real(pReal) Temperature,CPFEM_dt,J_inverse
 integer(pInt) CPFEM_mode               ! 1: regular computation, 2: collection, 3: recycling
 logical       CPFEM_updateJaco
!
 if (.not. CPFEM_init_done) then        ! initialization step
                                        ! three dimensional stress state check missing?
   call math_init()
   call mesh_init()
   call crystal_init()
   call constitutive_init()
   call CPFEM_init(Temperature)			
   CPFEM_init_done = .true.
 endif

 cp_en = mesh_FEasCP('elem',CPFEM_en)
  if (cp_en == 1 .and. CPFEM_in == 1) &
    write(6,'(a6,x,i4,x,a4,x,i4,x,a10,x,i2,x,a10,x,i2,x,a10,x,i2)') &
    'elem',cp_en,'IP',CPFEM_in,'theInc',theInc,'theCycle',theCycle,'theLovl',theLovl,'mode',CPFEM_mode

 select case (CPFEM_mode)
    case (2,1)     ! regular computation (with aging of results)
       if (.not. CPFEM_calc_done) then                ! puuh, me needs doing all the work...
           write (6,*) 'puuh me needs doing all the work', cp_en
           if (CPFEM_mode == 1) then                  ! age results at start of new increment
             CPFEM_Fp_old            = CPFEM_Fp_new
             constitutive_state_old  = constitutive_state_new
             write (6,*) '#### aged results'
           endif

           debug_cutbackDistribution = 0_pInt         ! initialize debugging data
           debug_InnerLoopDistribution = 0_pInt
           debug_OuterLoopDistribution = 0_pInt

           do e=1,mesh_NcpElems                       ! ## this shall be done in a parallel loop in the future ##
               do i=1,FE_Nips(mesh_element(2,e))      ! iterate over all IPs of this element's type
                  debugger = (e==1 .and. i==1)        ! switch on debugging for first IP in first element
                  call CPFEM_MaterialPoint(CPFEM_updateJaco, CPFEM_dt, i, e)
               enddo
           enddo

           call debug_info()                          ! output of debugging/performance statistics
           CPFEM_calc_done = .true.                   ! now calc is done
	   endif    
	   ! translate from P and dP/dF to CS and dCS/dE
       Kirchhoff_bar = matmul(CPFEM_PK1_bar(:,:,CPFEM_in, cp_en),transpose(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en)))
       J_inverse  = 1.0_pReal/math_det3x3(CPFEM_ffn1_bar(:,:,CPFEM_in, cp_en))
       CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel33to6(J_inverse*Kirchhoff_bar)
!
       H_bar = 0.0_pReal
       forall(i=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
         H_bar(i,j,k,l) = H_bar(i,j,k,l) + &
                          (CPFEM_ffn1_bar(j,m,CPFEM_in, cp_en)*CPFEM_ffn1_bar(l,n,CPFEM_in, cp_en)*CPFEM_dPdF_bar(i,m,k,n,CPFEM_in, cp_en) - &
                           math_I3(j,l)*CPFEM_ffn1_bar(i,m,CPFEM_in, cp_en)*CPFEM_PK1_bar(k,m, CPFEM_in, cp_en)) + &
                          0.5_pReal*(math_I3(i,k)*Kirchhoff_bar(j,l) + math_I3(j,l)*Kirchhoff_bar(i,k) + &
                                     math_I3(i,l)*Kirchhoff_bar(j,k) + math_I3(j,k)*Kirchhoff_bar(i,l))
       CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = math_Mandel3333to66(J_inverse*H_bar)

    case (3)    ! collect and return odd result
       CPFEM_Temperature(CPFEM_in,cp_en)  = Temperature
       CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)  = ffn
       CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) = ffn1
       CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_stress
       CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en) = CPFEM_odd_jacobian*math_identity2nd(CPFEM_ngens)
       CPFEM_calc_done = .false.

    case (4)    ! do nothing since we can recycle the former results (MARC specialty)
 end select
!
! return the local stress and the jacobian from storage
 CPFEM_stress(1:CPFEM_ngens) = CPFEM_stress_bar(1:CPFEM_ngens,CPFEM_in,cp_en)
 CPFEM_jaco(1:CPFEM_ngens,1:CPFEM_ngens) = CPFEM_jaco_bar(1:CPFEM_ngens,1:CPFEM_ngens,CPFEM_in,cp_en)
 if (cp_en == 1 .and. CPFEM_in == 1) write (6,*) 'stress',CPFEM_stress
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
 use prec, only: pReal,pInt,ijaco,nCutback
 use debug
 use math, only: math_pDecomposition,math_RtoEuler, inDeg, math_I3, math_invert3x3
 use IO,   only: IO_error
 use mesh, only: mesh_element
 use constitutive
 implicit none
!
 integer(pInt), parameter :: i_now = 1_pInt,i_then = 2_pInt
 character(len=128) msg
 integer(pInt) cp_en,CPFEM_in,grain,i,max_cutbacks
 logical updateJaco,error,cutback,post_flag
 real(pReal) CPFEM_dt,dt,t,volfrac,det
 real(pReal), dimension(3,3)     :: PK1
 real(pReal), dimension(3,3,3,3) :: dPdF
 real(pReal), dimension(3,3)     :: Fe,U,R,deltaFg,invFgthen,invFpnow,Lp
 real(pReal), dimension(3,3,2)   :: Fg,Fp
 real(pReal), dimension(constitutive_maxNstatevars,2) :: state
 real(pReal), dimension (:), allocatable              :: post_results
!
 CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) = 0.0_pReal                       ! zero out average first PK stress
 if (updateJaco) CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) = 0.0_pReal  ! zero out average consistent tangent
!
! -------------- grain loop -----------------
 do grain = 1,texture_Ngrains(mesh_element(4,cp_en))
   allocate(post_results(constitutive_Nresults(grain,CPFEM_in,cp_en))) ; post_results = 0.0_pReal
!
   i = 0_pInt                         ! cutback counter
   max_cutbacks = 0_pInt              ! maximum depth of cut backing
   dt = CPFEM_dt
   state(:,i_now) = constitutive_state_old(:,grain,CPFEM_in,cp_en)
   Fg(:,:,i_now)  = CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)
   Fp(:,:,i_now)  = CPFEM_Fp_old(:,:,grain,CPFEM_in,cp_en)
   invFgthen = 0.0_pReal
   invFpnow = 0.0_pReal
   call math_invert3x3(CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en),invFgthen,det,error)
   call math_invert3x3(Fp(:,:,i_now),invFpnow,det,error)
   if (dt /= 0.0_pReal) then
     Lp = (math_I3-matmul(Fp(:,:,i_now),matmul(invFgthen,matmul(Fg(:,:,i_now),invFpnow))))/dt ! fully plastic initial guess
   else
     Lp = 0.0_pReal                                          ! fully elastic guess 
   endif
!
   deltaFg = CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en)-CPFEM_ffn_bar(:,:,CPFEM_in,cp_en)
   Fg(:,:,i_then) = Fg(:,:,i_now)
   Fp(:,:,i_then) = Fp(:,:,i_now)
   state(:,i_then) = 0.0_pReal                               ! state_old as initial guess
   t = 0.0_pReal
   cutback = .false.                                         ! no cutback has happened so far
   msg = ''
   if (debugger) then
     write(6,*) 'required Fg from FEM'
     write(6,'(3(3(f5.3,x),/))') CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en)
     write(6,*) 'my Fp_old'
     write(6,'(3(3(f5.3,x),/))') CPFEM_Fp_old(:,:,grain,CPFEM_in,cp_en)
     write(6,*) 'my Fp_new'
     write(6,'(3(3(f5.3,x),/))') CPFEM_Fp_new(:,:,grain,CPFEM_in,cp_en)
     write(6,*) 'my state old'
     write(6,*) constitutive_state_old(:,grain,CPFEM_in,cp_en)
   endif

!
! ------- crystallite integration -----------
   do while ((t < CPFEM_dt) .or. (msg /= 'ok'))
!
     if (t+dt < CPFEM_dt) then          ! intermediate solution
       t = t+dt                         ! next time inc
       Fg(:,:,i_then) = Fg(:,:,i_then)+deltaFg  ! corresponding Fg
       post_flag = .false.
     else                               ! full step solution
       t = CPFEM_dt                     ! final time
       Fg(:,:,i_then) = CPFEM_ffn1_bar(:,:,CPFEM_in,cp_en) ! final Fg 
       post_flag = .true.
     endif
!
     if (debugger .and. CPFEM_dt > 0.0_pReal) &
       write (6,'(a,x,f7.5,x,a,x,f7.5,x,a,i2)') 'calculating from',(t-dt)/CPFEM_dt,'to',t/CPFEM_dt,'for grain',grain

     call CPFEM_Crystallite(msg,PK1,dPdF,post_results,post_flag,Lp,Fp(:,:,i_then),Fe,state(:,i_then),&
                            t,cp_en,CPFEM_in,grain,updateJaco .and. t==CPFEM_dt,&
                            Fg(:,:,i_then),Fp(:,:,i_now),state(:,i_now))
     if (msg == 'ok') then             ! solution converged
!       if (t == CPFEM_dt) then
!          debug_cutbackDistribution(max_cutbacks+1) = debug_cutbackDistribution(max_cutbacks+1)+1
!          exit                         ! reached final "then"
!       endif
       if (.not. cutback) then           ! stable solution at current speed?
         dt = 2.0_pReal*dt             ! double time-step
         i = i-1_pInt                  ! dec cutback counter
       endif
       cutback = .false.               ! solution in next step does not derive from a cutback
     else                              ! solution not found
       if (debugger) write (6,*) msg
       i = i+1_pInt                    ! inc cutback counter
       max_cutbacks = max(i,max_cutbacks)
       cutback = .true.
       if (i > nCutback) then          ! limit exceeded?
         debug_cutbackDistribution(nCutback+1) = debug_cutbackDistribution(nCutback+1)+1
         write(6,'(x,a,x,f10.8,x,a,x,f10.8,x,a,x,i6,x,a,x,i2,x,a,x,i2)') &
                 'inc fraction:',t/CPFEM_dt,'from',(t-dt)/CPFEM_dt,'element:',cp_en,'IP:',CPFEM_in,'grain:',grain
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

   debug_cutbackDistribution(max_cutbacks+1) = debug_cutbackDistribution(max_cutbacks+1)+1
!
! update crystallite matrices at t = t1
   CPFEM_Fp_new(:,:,grain,CPFEM_in,cp_en)         = Fp(:,:,i_then)
   constitutive_state_new(:,grain,CPFEM_in,cp_en) = state(:,i_then)
!
! contribute to IP result
    volfrac = constitutive_matVolFrac(grain,CPFEM_in,cp_en)*constitutive_texVolFrac(grain,CPFEM_in,cp_en)
    CPFEM_PK1_bar(:,:,CPFEM_in,cp_en) = CPFEM_PK1_bar(:,:,CPFEM_in,cp_en)+volfrac*PK1                  ! average Cauchy stress
    if (updateJaco) CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en) = CPFEM_dPdF_bar(:,:,:,:,CPFEM_in,cp_en)+volfrac*dPdF  ! consistent tangent
!
! update results plotted in MENTAT
   call math_pDecomposition(Fe,U,R,error) ! polar decomposition
   if (error) then
     write(6,*) Fe
     write(6,*) 'polar decomposition'
     write(6,*) 'Grain:             ',grain
     write(6,*) 'Integration point: ',CPFEM_in
     write(6,*) 'Element:           ',mesh_element(1,cp_en)
     call IO_error(650)
     return
   endif
   CPFEM_results(1:3,grain,CPFEM_in,cp_en) = math_RtoEuler(transpose(R))*inDeg        ! orientation
   CPFEM_results(4  ,grain,CPFEM_in,cp_en) = volfrac                                  ! volume fraction of orientation
   CPFEM_results(5:4+constitutive_Nresults(grain,CPFEM_in,cp_en),grain,CPFEM_in,cp_en) = post_results
!
   deallocate(post_results)
 enddo    ! grain loop
!
 return
!
 END SUBROUTINE
!
!********************************************************************
! Calculates the stress for a single component
!********************************************************************
 subroutine CPFEM_Crystallite(&
     msg,&        ! return message
     P,&          ! first PK stress
     dPdF,&       ! consistent tangent
     post_results,& ! plot results from constitutive model
     post_flag,&    ! its flag
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
!
 use prec, only: pReal,pInt,pert_Fg
 use debug
 use constitutive, only: constitutive_Nstatevars,constitutive_Nresults
 use mesh, only: mesh_element
 use math, only: math_Mandel6to33,math_Mandel33to6,math_Mandel3333to66,&
                 math_I3,math_det3x3,math_invert3x3
 implicit none
!
 character(len=*) msg
 logical updateJaco,error,post_flag
 integer(pInt) cp_en,CPFEM_in,grain,i,j,k,l,m,n
 real(pReal) dt,invJ,det
 real(pReal), dimension(3,3)     :: Lp,Lp_pert,Fg_old,Fg_new,Fg_pert,Fp_old,Fp_new,invFp_new,Fp_pert,invFp_pert
 real(pReal), dimension(3,3)     :: Fe_new,Fe_pert,Tstar,tau,P,P_pert,E_pert
 real(pReal), dimension(3,3,3,3) :: dPdF
 real(pReal), dimension(constitutive_Nstatevars(grain,CPFEM_in,cp_en)) :: state_old,state_new,state_pert
 real(pReal), dimension(constitutive_Nresults(grain,CPFEM_in,cp_en))   :: post_results
!
 call CPFEM_timeIntegration(msg,Lp,Fp_new,Fe_new,P,state_new,post_results,post_flag, &   ! def gradients and PK2 at end of time step
                            dt,cp_en,CPFEM_in,grain,Fg_new,Fp_old,state_old)
 if (msg /= 'ok') return                    ! solution not reached --> report back
 if (updateJaco) then                       ! consistent tangent using
                                            ! numerical perturbation of Fg (D. Tjahjanto Diss p.106)
   call math_invert3x3(Fp_new,invFp_new,det,error)
   if (error) then
     msg = 'inversion of Fp_new'
     return
   endif
   do k=1,3
     do l=1,3
       Fg_pert = Fg_new                       ! initialize perturbed Fg
       Fg_pert(k,l) = Fg_pert(k,l) + pert_Fg  ! perturb single component
       Lp_pert    = Lp
       state_pert = state_new                 ! initial guess from end of time step
       call CPFEM_timeIntegration(msg,Lp_pert,Fp_pert,Fe_pert,P_pert,state_pert,post_results,.false., &
                                  dt,cp_en,CPFEM_in,grain,Fg_pert,Fp_old,state_old)
       if (msg /= 'ok') then
         msg = 'consistent tangent --> '//msg
         return
       endif
!
       call math_invert3x3(Fp_pert,invFp_pert,det,error)
       if (error) then
         msg = 'inversion of Fp_pert'
         return
       endif
!
       dPdF(:,:,k,l) = (P_pert-P)/pert_Fg        ! constructin the tangent dP_ij/dFg_kl
     enddo
   enddo
 endif
!
 return
!
 END SUBROUTINE
!
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
     P,&                ! 1nd PK stress (taken as initial guess if /= 0)
     state,&            ! current microstructure at end of time inc (taken as guess if /= 0)
     results,&           ! post results from constitutive
     wantsConstitutiveResults,&        ! its flag
!
     dt,&               ! time increment
     cp_en,&            ! element number
     CPFEM_in,&         ! integration point number
     grain,&            ! grain number
     Fg_new,&           ! new total def gradient
     Fp_old,&           ! former plastic def gradient
     state_old)         ! former microstructure
!
 use prec
 use debug
 use mesh, only: mesh_element
 use constitutive, only: constitutive_Nstatevars,&
                         constitutive_homogenizedC,constitutive_dotState,constitutive_LpAndItsTangent,&
                         constitutive_Nresults,constitutive_Microstructure,constitutive_post_results
 use math
 implicit none
!
 character(len=*) msg
 logical failed,wantsConstitutiveResults
 integer(pInt) cp_en, CPFEM_in, grain
 integer(pInt) iOuter,iInner,dummy, i,j,k,l,m,n
 real(pReal) dt, det, p_hydro, leapfrog,maxleap
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(9,9) :: dLp,dTdLp,dRdLp,invdRdLp,eye2
 real(pReal), dimension(6,6) :: C_66
 real(pReal), dimension(3,3) :: Fg_new,invFg_new,Fp_new,invFp_new,Fp_old,invFp_old,Fe_new,Fe_old
 real(pReal), dimension(3,3) :: P,Tstar
 real(pReal), dimension(3,3) :: Lp,Lpguess,Lpguess_old,Rinner,Rinner_old,A,B,BT,AB,BTA
 real(pReal), dimension(3,3,3,3) :: C
 real(pReal), dimension(constitutive_Nstatevars(grain, CPFEM_in, cp_en)) :: state_old,state,ROuter
 real(pReal), dimension(constitutive_Nresults(grain,CPFEM_in,cp_en))   :: results
!
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
!
 Fe_old = matmul(Fg_new,invFp_old)
 A = matmul(transpose(Fe_old), Fe_old)
!
 if (all(state == 0.0_pReal)) state = state_old    ! former state guessed, if none specified
 iOuter = 0_pInt                                   ! outer counter
!
!
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
!
         iInner = 0_pInt
         leapfrog = 1.0_pReal                ! correction as suggested by invdRdLp-step
         maxleap = 1024.0_pReal              ! preassign maximum acceleration level
!
Inner: do              ! inner iteration: Lp
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
         if ((maxval(abs(Rinner)) < abstol_Inner) .or. &
             (any(abs(dt*Lpguess) > relevantStrain) .and. &
              maxval(abs(Rinner/Lpguess),abs(dt*Lpguess) > relevantStrain) < reltol_Inner)) &
           exit Inner
!
!          check for acceleration/deceleration in Newton--Raphson correction
         if (leapfrog > 1.0_pReal .and. &
             (sum(Rinner*Rinner) > sum(Rinner_old*Rinner_old) .or. &  ! worse residuum
              sum(Rinner*Rinner_old) < 0.0_pReal)) then               ! residuum changed sign (overshoot)
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
               if (debugger) then
                 write (6,*) msg
                 write (6,*) 'dRdLp',dRdLp
                 write (6,*) 'state',state
                 write (6,*) 'Lpguess',Lpguess
                 write (6,*) 'Tstar',Tstar_v
               endif                 
             return
           endif
!
           Rinner_old = Rinner                                        ! remember current residuum
           Lpguess_old = Lpguess                                      ! remember current Lp guess
           if (iInner > 1 .and. leapfrog < maxleap) leapfrog = 2.0_pReal * leapfrog   ! accelerate
         endif
!
         Lpguess = Lpguess_old                                        ! start from current guess
         Rinner  = Rinner_old                                         ! use current residuum
         forall (i=1:3,j=1:3,k=1:3,l=1:3) &                           ! leapfrog to updated Lpguess 
           Lpguess(i,j) = Lpguess(i,j) - leapfrog*invdRdLp(3*(i-1)+j,3*(k-1)+l)*Rinner(k,l)
       enddo Inner
!
       debug_InnerLoopDistribution(iInner) = debug_InnerLoopDistribution(iInner)+1
       ROuter = state - state_old - &
                dt*constitutive_dotState(Tstar_v,state,CPFEM_Temperature(CPFEM_in,cp_en),&
                                         grain,CPFEM_in,cp_en)          ! residuum from evolution of microstructure
       state = state - ROuter                                           ! update of microstructure
       if (maxval(abs(Router/state),state /= 0.0_pReal) < reltol_Outer) exit Outer
     enddo Outer
!
 debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
 invFp_new = matmul(invFp_old,B)
 call math_invert3x3(invFp_new,Fp_new,det,failed)
 if (failed) then
   msg = 'inversion Fp_new'
   return
 endif
!
 if (wantsConstitutiveResults) then     ! get the post_results upon request
   results = 0.0_pReal
   results = constitutive_post_results(Tstar_v,state,dt,CPFEM_Temperature(CPFEM_in,cp_en),grain,CPFEM_in,cp_en)
 endif
!
 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal)     ! regularize Fp by det = det(InvFp_new) !!
 Fe_new = matmul(Fg_new,invFp_new)              ! calc resulting Fe
 forall (i=1:3) Tstar_v(i) = Tstar_v(i)+p_hydro ! add hydrostatic component back
 P = matmul(Fe_new,matmul(Tstar,transpose(invFp_new)))    ! first PK stress
!
 return
!
 END SUBROUTINE
!
 END MODULE
!##############################################################

