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
     Fp_new,&     ! new plastic deformation gradient
     Fe_new,&     ! new "elastic" deformation gradient
     state_new,&  ! new state variable array
!
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
 use prec, only: pReal,pInt,pert_Fg,subStepMin
 use debug
 use constitutive, only: constitutive_Nstatevars,constitutive_Nresults
 use mesh, only: mesh_element
 use math, only: math_Mandel6to33,math_Mandel33to6,math_Mandel3333to66,&
                 math_I3,math_det3x3,math_invert3x3
 implicit none
!
 character(len=*) msg
 logical updateJaco,error,cuttedBack,guessNew
 integer(pInt) cp_en,ip,grain,i,j,k,l,m,n, nCutbacks
 real(pReal) Temperature
 real(pReal) dt,dt_aim,subFrac,subStep,invJ,det
 real(pReal), dimension(3,3)     :: Lp,Lp_pert,inv
 real(pReal), dimension(3,3)     :: Fg_old,Fg_current,Fg_aim,Fg_new,Fg_pert,deltaFg
 real(pReal), dimension(3,3)     :: Fp_old,Fp_current,Fp_new,Fp_pert
 real(pReal), dimension(3,3)     :: Fe_old,Fe_current,Fe_new,Fe_pert
 real(pReal), dimension(3,3)     :: Tstar,tau,P,P_pert
 real(pReal), dimension(3,3,3,3) :: dPdF
 real(pReal), dimension(constitutive_Nstatevars(grain,ip,cp_en)) :: state_old,state_current,state_new,state_pert
 real(pReal), dimension(constitutive_Nresults(grain,ip,cp_en))   :: post_results
!
 deltaFg = Fg_new-Fg_old
 subFrac = 0.0_pReal
 subStep = 1.0_pReal
 nCutbacks = 0_pInt
!
 Fg_aim = Fg_old                                                   ! make "new", "aim" a synonym for "old"
 Fp_new = Fp_old
 call math_invert3x3(Fp_old,inv,det,error)
 Fe_new = matmul(Fg_old,inv)
 state_new = state_old
!
 cuttedBack = .false.
 guessNew = .true.
!
! begin the cutback loop
 do while (subStep > subStepMin)                                   ! continue until finished or too much cut backing
   if (.not. cuttedBack) then
     Fg_current = Fg_aim                                           ! wind forward
     Fp_current = Fp_new
     Fe_current = Fe_new
     state_current = state_new
   endif
!
   Fg_aim = Fg_current + subStep*deltaFg                           ! aim for Fg
   dt_aim = subStep*dt                                             ! aim for dt
   msg = ''                                                        ! error free so far
   if (guessNew) then                                              ! calculate new Lp guess when cutted back
     if (dt_aim /= 0.0_pReal) then
       call math_invert3x3(Fg_aim,inv,det,error)
       Lp = (math_I3-matmul(Fp_current,matmul(inv,Fe_current)))/dt ! fully plastic initial guess
     else
       Lp = 0.0_pReal                                              ! fully elastic guess 
     endif
   endif
   call TimeIntegration(msg,Lp,Fp_new,Fe_new,P,state_new,post_results,.true., &   ! def gradients and PK2 at end of time step
                        dt_aim,cp_en,ip,grain,Temperature,Fg_aim,Fp_current,state_current)
!
   if (msg == 'ok') then
     cuttedBack = .false.                 ! no cut back required
     guessNew = .false.                   ! keep the Lp
     subFrac = subFrac + subStep
     subStep = 1.0_pReal - subFrac        ! try one go
   else
     nCutbacks = nCutbacks + 1            ! record additional cutback
     cuttedBack = .true.                  ! encountered problems -->
     guessNew = .true.                    ! redo plastic Lp guess
     subStep = subStep / 2.0_pReal        ! cut time step in half
   endif
 enddo
!
 debug_cutbackDistribution(nCutbacks+1) = debug_cutbackDistribution(nCutbacks+1)+1
!
 if (msg /= 'ok') return                    ! solution not reached --> report back
 if (updateJaco) then                       ! consistent tangent using
   do k=1,3
     do l=1,3
       Fg_pert = Fg_new                       ! initialize perturbed Fg
       Fg_pert(k,l) = Fg_pert(k,l) + pert_Fg  ! perturb single component
       Lp_pert    = Lp
       state_pert = state_new                 ! initial guess from end of time step
       call TimeIntegration(msg,Lp,Fp_pert,Fe_pert,P_pert,state_pert,post_results,.false., &   ! def gradients and PK2 at end of time step
                                  dt_aim,cp_en,ip,grain,Temperature,Fg_pert,Fp_current,state_current)
       if (msg /= 'ok') then
         msg = 'consistent tangent --> '//msg
         return
       endif
       dPdF(:,:,k,l) = (P_pert-P)/pert_Fg        ! constructing the tangent dP_ij/dFg_kl from forward differences
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
 integer(pInt) cp_en, ip, grain
 integer(pInt) iOuter,iInner,dummy, i,j,k,l,m,n
 real(pReal) dt, Temperature, det, p_hydro, leapfrog,maxleap
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(9,9) :: dLp,dTdLp,dRdLp,invdRdLp,eye2
 real(pReal), dimension(6,6) :: C_66
 real(pReal), dimension(3,3) :: Fg_new,Fp_new,invFp_new,Fp_old,invFp_old,Fe_new,Fe_old
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

 A = matmul(transpose(invFp_old), matmul(transpose(Fg_new),matmul(Fg_new,invFp_old)))
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
         call constitutive_Microstructure(state,Temperature,grain,ip,cp_en)
         C_66 = constitutive_HomogenizedC(state, grain, ip, cp_en)
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
                                           Tstar_v,state,Temperature,grain,ip,cp_en)
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
                 write (6,'(a,/,9(9(e9.3,x)/))') 'dRdLp', dRdLp(1:9,:)
                 write (6,*) 'state',state
                 write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
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
                dt*constitutive_dotState(Tstar_v,state,Temperature,&
                                         grain,ip,cp_en)          ! residuum from evolution of microstructure
       state = state - ROuter                                           ! update of microstructure
       if (maxval(abs(Router/state),state /= 0.0_pReal) < reltol_Outer) exit Outer
     enddo Outer
!
 debug_OuterLoopDistribution(iOuter) = debug_OuterLoopDistribution(iOuter)+1


 invFp_new = matmul(invFp_old,B)
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
 Fe_new = matmul(Fg_new,invFp_new)              ! calc resulting Fe
 forall (i=1:3) Tstar_v(i) = Tstar_v(i)+p_hydro ! add hydrostatic component back
 P = matmul(Fe_new,matmul(Tstar,transpose(invFp_new)))    ! first PK stress
!
 return
!
 END SUBROUTINE
!!
 END MODULE
!##############################################################

