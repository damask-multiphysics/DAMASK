!##############################################################
 MODULE crystallite
!##############################################################
!    *** Solution at single crystallite level ***
!
CONTAINS
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
 subroutine SingleCrystallite(&
     msg,&        ! return message
     P,&          ! first PK stress
     dPdF,&       ! consistent tangent
     post_results,& ! plot results from constitutive model
     Lp,&         ! plastic velocity gradient
     Fp_new,&     ! new plastic deformation gradient
     Fe_new,&     ! new "elastic" deformation gradient
     state_new,&  ! new state variable array
     dt,&         ! time increment
     cp_en,&      ! element number
     ip,&         ! integration point number
     grain,&      ! grain number
     updateJaco,& ! update of Jacobian required
     Temperature,& ! temperature of crystallite
     Fg_new,&     ! new global deformation gradient
     Fg_old,&     ! old global deformation gradient
     Fp_old,&     ! old plastic deformation gradient
     state_old)   ! old state variable array
!
 use prec, only: pReal,pInt,pert_Fg,subStepMin, nCutback
 use debug
 use constitutive, only: constitutive_Nstatevars,constitutive_Nresults
 use mesh, only: mesh_element
 use math

 implicit none
!
 character(len=*) msg
 logical updateJaco,error,success,guessNew
 integer(pInt) cp_en,ip,grain,k,l, nCutbacks, maxCutbacks
 real(pReal) Temperature
 real(pReal) dt,dt_aim,subFrac,subStep,det
 real(pReal), dimension(3,3)     :: Lp,Lp_interpolated,Lp_pert,inv
 real(pReal), dimension(3,3)     :: Fg_old,Fg_current,Fg_new,Fg_pert,Fg_aim,deltaFg
 real(pReal), dimension(3,3)     :: Fp_old,Fp_current,Fp_new,Fp_pert
 real(pReal), dimension(3,3)     :: Fe_current,Fe_new,Fe_pert
 real(pReal), dimension(3,3)     :: P,P_pert
 real(pReal), dimension(3,3,3,3) :: dPdF
 real(pReal), dimension(constitutive_Nstatevars(grain,ip,cp_en)) :: state_old,state_new
 real(pReal), dimension(constitutive_Nstatevars(grain,ip,cp_en)) :: state_current,state_bestguess,state_pert
 real(pReal), dimension(constitutive_Nresults(grain,ip,cp_en))   :: post_results
 
 deltaFg = Fg_new - Fg_old
 subFrac = 0.0_pReal
 subStep = 1.0_pReal
 nCutbacks = 0_pInt
 maxCutbacks = 0_pInt
 Fg_current = Fg_old                                               ! initialize to start of inc
 Fp_current = Fp_old
 call math_invert3x3(Fp_old,inv,det,error)
 Fe_current = math_mul33x33(Fg_old,inv)
 state_current = state_old
 success = .false.                                                 ! pretend cutback
 dt_aim = 0.0_pReal                                                ! prevent initial Lp interpolation
 
! begin the cutback loop
 do while (subStep > subStepMin)                                   ! continue until finished or too much cut backing
   if (success) then                                               ! wind forward
     Fg_current = Fg_aim
	 Fe_current = Fe_new
	 Fp_current = Fp_new
	 state_current = state_new
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

   Fg_aim = Fg_current + subStep*deltaFg                           ! aim for Fg
   dt_aim = subStep*dt                                             ! aim for dt
   if (debugger) then
!$OMP CRITICAL (write2out)
     write (6,*) 'using these values'
     write (6,'(a,/,3(4(f9.3,x)/))') 'state current / MPa',state_current/1e6_pReal
     write (6,'(a,/,3(4(f9.3,x)/))') 'state new / MPa',state_new/1e6_pReal
     write (6,'(a,/,3(3(f12.7,x)/))') 'Fe current',Fe_current(1:3,:)
     write (6,'(a,/,3(3(f12.7,x)/))') 'Fp current',Fp_current(1:3,:)
     write (6,'(a,/,3(3(f12.7,x)/))') 'Lp (old=new guess)',Lp(1:3,:)
     write (6,'(a20,f,x,a2,x,f)') 'integrating from ',subFrac,'to',(subFrac+subStep)
!$OMP END CRITICAL (write2out)
   endif
   call TimeIntegration(msg,Lp,Fp_new,Fe_new,P,state_new,post_results, &   ! def gradients and PK2 at end of time step
                        maxval(abs(Fg_aim-Fg_new)) < relevantStrain, &     ! post results only if asking for final values
                        dt_aim,cp_en,ip,grain,Temperature,Fg_aim,Fp_current,state_current)

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
	 if (debugger) then
!$OMP CRITICAL (write2out)
	   write (6,*) '>>>>>>>>>>>>>>>>>>>> cutback <<<<<<<<<<<<<<<<<<<<<<'
!$OMP END CRITICAL (write2out)
     endif

   endif
 enddo  ! potential substepping
!
!$OMP CRITICAL (cutback)
 debug_cutbackDistribution(min(nCutback,maxCutbacks)+1) = debug_cutbackDistribution(min(nCutback,maxCutbacks)+1)+1
!$OMP END CRITICAL (cutback)
!
 if (msg /= 'ok') return                    ! solution not reached --> report back
 if (updateJaco) then                       ! consistent tangent using
   if (debugger) then
!$OMP CRITICAL (write2out)
     write (6,*) 'Jacobian calc'
!$OMP END CRITICAL (write2out)
   endif   
   do k=1,3
     do l=1,3
       Fg_pert = Fg_new                       ! initialize perturbed Fg
       Fg_pert(k,l) = Fg_pert(k,l) + pert_Fg  ! perturb single component
       Lp_pert    = Lp
       state_pert = state_new                 ! initial guess from end of time step
       call TimeIntegration(msg,Lp_pert,Fp_pert,Fe_pert,P_pert,state_pert,post_results,.false., &   ! def gradients and PK2 at end of time step
                            dt_aim,cp_en,ip,grain,Temperature,Fg_pert,Fp_current,state_current)

       if (msg == 'ok') &
         dPdF(:,:,k,l) = (P_pert-P)/pert_Fg   ! constructing tangent dP_ij/dFg_kl only if valid forward difference
                                              ! otherwise leave component unchanged
     enddo
   enddo
 endif
!

 msg = 'ok' ! a new consistent tangent was computed even if msg was not ok for all components

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
 SUBROUTINE TimeIntegration(&
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
     ip,&               ! integration point number
     grain,&            ! grain number
     Temperature,&      ! temperature
     Fg_new,&           ! new total def gradient
     Fp_old,&           ! former plastic def gradient
     state_old)         ! former microstructure
 use prec
 use debug
 use mesh, only: mesh_element
 use constitutive, only: constitutive_Nstatevars,&
                         constitutive_homogenizedC,constitutive_dotState,constitutive_LpAndItsTangent,&
                         constitutive_Nresults,constitutive_Microstructure,constitutive_post_results
 use math
 use IO
 implicit none
!
 character(len=*) msg
 logical failed,wantsConstitutiveResults
 integer(pInt) cp_en, ip, grain
 integer(pInt) iOuter,iInner,dummy, i,j,k,l,m,n
 real(pReal) dt, Temperature, det, p_hydro, leapfrog,maxleap
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(9,9) :: dLp,dTdLp,dRdLp,invdRdLp,eye2
 real(pReal), dimension(6,6) :: C_66
 real(pReal), dimension(3,3) :: Fg_new,Fp_new,invFp_new,Fp_old,invFp_old,Fe_new
 real(pReal), dimension(3,3) :: P,Tstar
 real(pReal), dimension(3,3) :: Lp,Lpguess,Lpguess_old,Rinner,Rinner_old,A,B,BT,AB,BTA
 real(pReal), dimension(3,3,3,3) :: C
 real(pReal), dimension(constitutive_Nstatevars(grain, ip, cp_en)) :: state_old,state,ROuter
 real(pReal), dimension(constitutive_Nresults(grain,ip,cp_en))   :: results
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
 if (all(state == 0.0_pReal)) state = state_old    ! former state guessed, if none specified
 iOuter = 0_pInt                                   ! outer counter
!
 if (debugger) then
!$OMP CRITICAL (write2out)
   write (6,'(a,/,3(3(f12.7,x)/))') 'Fg to be calculated',Fg_new
!$OMP END CRITICAL (write2out)
 endif
!
Outer: do                ! outer iteration: State
         iOuter = iOuter+1
         if (debugger) then
!$OMP CRITICAL (write2out)
           write (6,'(a,i3)') '---outer ',iOuter
           write (6,'(a,/,3(4(f9.3,x)/))') 'state old / MPa',state_old/1e6_pReal
           write (6,'(a,/,3(4(f9.3,x)/))') 'state / MPa',state/1e6_pReal
           write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
!$OMP END CRITICAL (write2out)
         endif

         if (iOuter > nOuter) then
           msg = 'limit Outer iteration'
!$OMP CRITICAL (out)
           debug_OuterLoopDistribution(nOuter) = debug_OuterLoopDistribution(nOuter)+1
!$OMP END CRITICAL (out)
           return
         endif
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
         if (debugger) then
!$OMP CRITICAL (write2out)
           write (6,'(a,i3)') 'inner ',iInner
		   if (wantsConstitutiveResults .and. iOuter == 1 .and. iInner < 3) then
             write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
		   endif
!$OMP END CRITICAL (write2out)
         endif
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
         Tstar = math_Mandel6to33(Tstar_v)
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
       ROuter = state - state_old - &
                dt*constitutive_dotState(Tstar_v,state,Temperature,&
                                         grain,ip,cp_en)          ! residuum from evolution of microstructure
       state = state - ROuter                                           ! update of microstructure

	   if (iOuter==nOuter) then
!!$OMP CRITICAL (write2out)
	      write (6,*) 'Terminated outer loop at el,ip,grain',cp_en,ip,grain
!!$OMP END CRITICAL (write2out)
	      exit Outer
	   endif
       if (maxval(abs(Router/state),state /= 0.0_pReal) < reltol_Outer) exit Outer
     enddo Outer
!
!$OMP CRITICAL (out)
 debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1
!$OMP END CRITICAL (out)
 invFp_new = math_mul33x33(invFp_old,B)
 call math_invert3x3(invFp_new,Fp_new,det,failed)
 if (failed) then
   msg = 'inversion Fp_new^-1'
   return
 endif
!
 if (wantsConstitutiveResults) then     ! get the post_results upon request
   results = 0.0_pReal
   results = constitutive_post_results(Tstar_v,state,Temperature,dt,grain,ip,cp_en)
 endif
!
 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal)     ! regularize Fp by det = det(InvFp_new) !!
 forall (i=1:3) Tstar_v(i) = Tstar_v(i)+p_hydro ! add hydrostatic component back
 Fe_new = math_mul33x33(Fg_new,invFp_new)              ! calc resulting Fe
 P = math_mul33x33(Fe_new,math_mul33x33(Tstar,transpose(invFp_new)))    ! first PK stress

 return
!
 END SUBROUTINE
!
!
 END MODULE
!##############################################################

