
!***************************************
!*      Module: CRYSTALLITE            *
!***************************************
!* contains:                           *
!* - _init                             *
!* - materialpoint_stressAndItsTangent *
!* - _partitionDeformation             *
!* - _updateState                      *
!* - _averageStressAndItsTangent       *
!* - _postResults                      *
!***************************************

MODULE crystallite

 use prec, only: pReal,pInt
 implicit none
!
! ****************************************************************
! *** General variables for the crystallite calculation        ***
! ****************************************************************
 integer(pInt), parameter :: crystallite_Nresults = 5_pInt                        ! phaseID, volfrac within this phase, Euler angles 

 real(pReal), dimension (:,:,:,:,:),    allocatable :: crystallite_Fe, &          ! current "elastic" def grad (end of converged time step)
                                                       crystallite_Fp, &          ! current plastic def grad (end of converged time step)
                                                       crystallite_Lp, &          ! current plastic velocitiy grad (end of converged time step)
                                                       crystallite_F0, &          ! def grad at start of FE inc
                                                       crystallite_Fp0, &         ! plastic def grad at start of FE inc
                                                       crystallite_Lp0, &         ! plastic velocitiy grad at start of FE inc
                                                       crystallite_partionedF,  & ! def grad to be reached at end of homog inc
                                                       crystallite_partionedF0, & ! def grad at start of homog inc
                                                       crystallite_partionedFp0,& ! plastic def grad at start of homog inc
                                                       crystallite_partionedLp0,& ! plastic velocity grad at start of homog inc
                                                       crystallite_subF,  &       ! def grad to be reached at end of crystallite inc
                                                       crystallite_subF0, &       ! def grad at start of crystallite inc
                                                       crystallite_subFp0,&       ! plastic def grad at start of crystallite inc
                                                       crystallite_subLp0,&       ! plastic velocity grad at start of crystallite inc
                                                       crystallite_P              ! 1st Piola-Kirchhoff stress per grain
 real(pReal), dimension (:,:,:,:),      allocatable :: crystallite_Tstar_v        ! 2nd Piola-Kirchhoff stress (vector) per grain
 real(pReal), dimension (:,:,:,:,:,:,:),allocatable :: crystallite_dPdF, &        ! individual dPdF per grain
                                                       crystallite_fallbackdPdF   ! dPdF fallback for non-converged grains (elastic prediction)
 real(pReal), dimension (:,:,:),        allocatable :: crystallite_dt, &          ! requested time increment of each grain
                                                       crystallite_subdt, &       ! substepped time increment of each grain
                                                       crystallite_subFrac, &     ! already calculated fraction of increment
                                                       crystallite_subStep, &     ! size of next integration step
                                                       crystallite_Temperature    ! Temp of each grain

 logical, dimension (:,:,:),      allocatable :: crystallite_localConstitution, & ! indicates this grain to have purely local constitutive law
                                                 crystallite_requested, &         ! flag to request crystallite calculation
                                                 crystallite_onTrack, &           ! flag to indicate ongoing calculation
                                                 crystallite_converged            ! convergence flag


 CONTAINS

!********************************************************************
! allocate and initialize per grain variables
!********************************************************************
 subroutine crystallite_init()

 use prec,      only: pInt,pReal
 use debug,     only: debug_info,debug_reset
 use math,      only: math_I3,math_EulerToR
 use FEsolving, only: FEsolving_execElem,FEsolving_execIP
 use mesh,      only: mesh_element,mesh_NcpElems,mesh_maxNips
 use material,  only: homogenization_Ngrains,homogenization_maxNgrains,&
                      material_EulerAngles,material_phase,phase_localConstitution
 implicit none

 integer(pInt) g,i,e, gMax,iMax,eMax, myNgrains
 
 gMax = homogenization_maxNgrains
 iMax = mesh_maxNips
 eMax = mesh_NcpElems
 
 allocate(crystallite_Fe(3,3,gMax,iMax,eMax));                         crystallite_Fe = 0.0_pReal
 allocate(crystallite_Fp(3,3,gMax,iMax,eMax));                         crystallite_Fp = 0.0_pReal
 allocate(crystallite_Lp(3,3,gMax,iMax,eMax));                         crystallite_Lp = 0.0_pReal
 allocate(crystallite_F0(3,3,gMax,iMax,eMax));                         crystallite_F0 = 0.0_pReal
 allocate(crystallite_Fp0(3,3,gMax,iMax,eMax));                       crystallite_Fp0 = 0.0_pReal
 allocate(crystallite_Lp0(3,3,gMax,iMax,eMax));                       crystallite_Lp0 = 0.0_pReal
 allocate(crystallite_partionedF(3,3,gMax,iMax,eMax));        crystallite_partionedF0 = 0.0_pReal
 allocate(crystallite_partionedF0(3,3,gMax,iMax,eMax));       crystallite_partionedF0 = 0.0_pReal
 allocate(crystallite_partionedFp0(3,3,gMax,iMax,eMax));     crystallite_partionedFp0 = 0.0_pReal
 allocate(crystallite_partionedLp0(3,3,gMax,iMax,eMax));     crystallite_partionedLp0 = 0.0_pReal
 allocate(crystallite_subF(3,3,gMax,iMax,eMax));                     crystallite_subF = 0.0_pReal
 allocate(crystallite_subF0(3,3,gMax,iMax,eMax));                   crystallite_subF0 = 0.0_pReal
 allocate(crystallite_subFp0(3,3,gMax,iMax,eMax));                 crystallite_subFp0 = 0.0_pReal
 allocate(crystallite_subLp0(3,3,gMax,iMax,eMax));                 crystallite_subLp0 = 0.0_pReal
 allocate(crystallite_P(3,3,gMax,iMax,eMax));                           crystallite_P = 0.0_pReal
 allocate(crystallite_Tstar_v(6,gMax,iMax,eMax));                 crystallite_Tstar_v = 0.0_pReal
 allocate(crystallite_dPdF(3,3,3,3,gMax,iMax,eMax));                 crystallite_dPdF = 0.0_pReal
 allocate(crystallite_fallbackdPdF(3,3,3,3,gMax,iMax,eMax)); crystallite_fallbackdPdF = 0.0_pReal
 allocate(crystallite_dt(gMax,iMax,eMax));                             crystallite_dt = 0.0_pReal
 allocate(crystallite_subdt(gMax,iMax,eMax));                       crystallite_subdt = 0.0_pReal
 allocate(crystallite_subFrac(gMax,iMax,eMax));                   crystallite_subFrac = 0.0_pReal
 allocate(crystallite_subStep(gMax,iMax,eMax));                   crystallite_subStep = 0.0_pReal
 allocate(crystallite_Temperature(gMax,iMax,eMax));           crystallite_Temperature = 0.0_pReal
 allocate(crystallite_localConstitution(gMax,iMax,eMax));
 allocate(crystallite_requested(gMax,iMax,eMax));               crystallite_requested = .false.
 allocate(crystallite_onTrack(gMax,iMax,eMax));                   crystallite_onTrack = .false.
 allocate(crystallite_converged(gMax,iMax,eMax));               crystallite_converged = .true.

!$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over all cp elements
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element
     do g = 1,myNgrains
       crystallite_Fp0(:,:,g,i,e) = math_EulerToR(material_EulerAngles(:,g,i,e))          ! plastic def gradient reflects init orientation
       crystallite_F0(:,:,g,i,e)  = math_I3
       crystallite_partionedFp0(:,:,g,i,e) = crystallite_Fp0(:,:,g,i,e)
       crystallite_partionedF0(:,:,g,i,e)  = crystallite_F0(:,:,g,i,e)
       crystallite_partionedF(:,:,g,i,e)   = crystallite_F0(:,:,g,i,e)
       crystallite_requested(g,i,e)        = .true.
       crystallite_localConstitution(g,i,e) = phase_localConstitution(material_phase(g,i,e))
     enddo
   enddo
 enddo
!$OMP END PARALLEL DO

 call crystallite_stressAndItsTangent(.true.)                 ! request elastic answers
 crystallite_fallbackdPdF = crystallite_dPdF                  ! use initial elastic stiffness as fallback
 
!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  crystallite init  -+>>>'
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Nresults:          ', crystallite_Nresults
 write(6,'(a32,x,7(i5,x))') 'crystallite_Fe:                ', shape(crystallite_Fe)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Fp:                ', shape(crystallite_Fp)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Lp:                ', shape(crystallite_Lp)
 write(6,'(a32,x,7(i5,x))') 'crystallite_F0:                ', shape(crystallite_F0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Fp0:               ', shape(crystallite_Fp0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Lp0:               ', shape(crystallite_Lp0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_partionedF:        ', shape(crystallite_partionedF)
 write(6,'(a32,x,7(i5,x))') 'crystallite_partionedF0:       ', shape(crystallite_partionedF0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_partionedFp0:      ', shape(crystallite_partionedFp0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_partionedLp0:      ', shape(crystallite_partionedLp0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subF:              ', shape(crystallite_subF)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subF0:             ', shape(crystallite_subF0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subFp0:            ', shape(crystallite_subFp0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subLp0:            ', shape(crystallite_subLp0)
 write(6,'(a32,x,7(i5,x))') 'crystallite_P:                 ', shape(crystallite_P)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Tstar_v:           ', shape(crystallite_Tstar_v)
 write(6,'(a32,x,7(i5,x))') 'crystallite_dPdF:              ', shape(crystallite_dPdF)
 write(6,'(a32,x,7(i5,x))') 'crystallite_fallbackdPdF:      ', shape(crystallite_fallbackdPdF)
 write(6,'(a32,x,7(i5,x))') 'crystallite_dt:                ', shape(crystallite_dt)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subdt:             ', shape(crystallite_subdt)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subFrac:           ', shape(crystallite_subFrac)
 write(6,'(a32,x,7(i5,x))') 'crystallite_subStep:           ', shape(crystallite_subStep)
 write(6,'(a32,x,7(i5,x))') 'crystallite_Temperature:       ', shape(crystallite_Temperature)
 write(6,'(a32,x,7(i5,x))') 'crystallite_localConstitution: ', shape(crystallite_localConstitution)
 write(6,'(a32,x,7(i5,x))') 'crystallite_requested:         ', shape(crystallite_requested)
 write(6,'(a32,x,7(i5,x))') 'crystallite_onTrack:           ', shape(crystallite_onTrack)
 write(6,'(a32,x,7(i5,x))') 'crystallite_converged:         ', shape(crystallite_converged)
 write(6,*)
 write(6,*) 'Number of non-local grains: ',count(.not. crystallite_localConstitution)
 call flush(6)
!$OMP END CRITICAL (write2out)

 call debug_info()
 call debug_reset()
 
 return

 end subroutine


!********************************************************************
! calculate stress (P) and tangent (dPdF) for crystallites
!********************************************************************
subroutine crystallite_stressAndItsTangent(updateJaco)

 use prec,         only: pInt,pReal,subStepMin,nCryst
 use debug
 use IO,           only: IO_warning
 use math
 use FEsolving,    only: FEsolving_execElem, FEsolving_execIP
 use mesh,         only: mesh_element
 use material,     only: homogenization_Ngrains
 use constitutive
 implicit none
 
 logical, intent(in) :: updateJaco
 real(pReal), dimension(3,3) :: invFp,Fe_guess,PK2,myF,myFp,myFe,myLp,myP
 real(pReal), dimension(constitutive_maxSizeState) :: myState
 integer(pInt) crystallite_Niteration
 integer(pInt) g,i,e,k,l, myNgrains, mySizeState
 logical converged

! ------ initialize to starting condition ------

 write (6,*)
 write (6,*) 'Crystallite request from Materialpoint'
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedF0  of 1 8 1',crystallite_partionedF0(1:3,:,1,8,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedFp0 of 1 8 1',crystallite_partionedFp0(1:3,:,1,8,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedF   of 1 8 1',crystallite_partionedF(1:3,:,1,8,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedLp0 of 1 8 1',crystallite_partionedLp0(1:3,:,1,8,1)


!$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
     do g = 1,myNgrains
       if (crystallite_requested(g,i,e)) then                 ! initialize restoration point of ...
         constitutive_subState0(g,i,e)%p = constitutive_partionedState0(g,i,e)%p  ! ...microstructure
         crystallite_subFp0(:,:,g,i,e)   = crystallite_partionedFp0(:,:,g,i,e)    ! ...plastic def grad
         crystallite_subLp0(:,:,g,i,e)   = crystallite_partionedLp0(:,:,g,i,e)    ! ...plastic velocity grad
         crystallite_subF0(:,:,g,i,e)    = crystallite_partionedF0(:,:,g,i,e)     ! ...def grad

         crystallite_subFrac(g,i,e) = 0.0_pReal
         crystallite_subStep(g,i,e) = 2.0_pReal
         crystallite_onTrack(g,i,e) = .true.
         crystallite_converged(g,i,e) = .false.               ! pretend failed step of twice the required size
       endif
     enddo
   enddo
 enddo
!$OMP END PARALLEL DO

! ------ cutback loop ------

 do while (any(crystallite_subStep(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMin))  ! cutback loop for crystallites
   write(6,*)
   write(6,*) 'entering cutback at crystallite_stress'
   if (any(.not. crystallite_converged(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) .and. &         ! any non-converged grain
           .not. crystallite_localConstitution(:,:,FEsolving_execELem(1):FEsolving_execElem(2))) ) &    ! has non-local constitution?
     crystallite_converged(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) = &
       crystallite_converged(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) .and. &
       crystallite_localConstitution(:,:,FEsolving_execELem(1):FEsolving_execElem(2))       ! reset non-local grains' convergence status
!$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
       do g = 1,myNgrains
         if (crystallite_converged(g,i,e)) then
           crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
           crystallite_subStep(g,i,e) = min(1.0_pReal-crystallite_subFrac(g,i,e), 2.0_pReal * crystallite_subStep(g,i,e))
           if (crystallite_subStep(g,i,e) > subStepMin) then
             crystallite_subF0(:,:,g,i,e)    = crystallite_subF(:,:,g,i,e)         ! wind forward...
             crystallite_subFp0(:,:,g,i,e)   = crystallite_Fp(:,:,g,i,e)           ! ...plastic def grad
             crystallite_subLp0(:,:,g,i,e)   = crystallite_Lp(:,:,g,i,e)           ! ...plastic velocity gradient
             constitutive_subState0(g,i,e)%p = constitutive_state(g,i,e)%p         ! ...microstructure
           endif
         else
           crystallite_subStep(g,i,e) = 0.5_pReal * crystallite_subStep(g,i,e)     ! cut step in half and restore...
           crystallite_Fp(:,:,g,i,e) = crystallite_subFp0(:,:,g,i,e)               ! ...plastic def grad
           crystallite_Lp(:,:,g,i,e) = crystallite_subLp0(:,:,g,i,e)               ! ...plastic velocity grad
           constitutive_state(g,i,e)%p = constitutive_subState0(g,i,e)%p           ! ...microstructure
         endif
     
         crystallite_onTrack(g,i,e) = crystallite_subStep(g,i,e) > subStepMin      ! still on track or already done (beyond repair)
         if (crystallite_onTrack(g,i,e)) then                                      ! specify task (according to substep)
           crystallite_subF(:,:,g,i,e)  = crystallite_subF0(:,:,g,i,e) + &
                                          crystallite_subStep(g,i,e) * &
                                          (crystallite_partionedF(:,:,g,i,e) - crystallite_partionedF0(:,:,g,i,e))
           crystallite_subdt(g,i,e)     = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
           crystallite_converged(g,i,e) = .false.                                  ! start out non-converged
         endif
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO

! ------ convergence loop for stress and state ------

   crystallite_Niteration = 0_pInt
   
   do while (any(            crystallite_requested(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and.       crystallite_onTrack(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. crystallite_converged(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. crystallite_Niteration < nCryst)              ! convergence loop for crystallite
     crystallite_Niteration = crystallite_Niteration + 1


! --+>> stress integration <<+--
!
! incrementing by crystallite_subdt
! based on crystallite_subF0,.._subFp0,.._subLp0
!          constitutive_state is internally interpolated with .._subState0
!          to account for substepping within _integrateStress
! results in crystallite_Fp,.._Lp

!$OMP PARALLEL DO
     do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
         do g = 1,myNgrains
           if (crystallite_requested(g,i,e) .and. &
               crystallite_onTrack(g,i,e)) &                      ! all undone crystallites
             crystallite_onTrack(g,i,e) = crystallite_integrateStress(g,i,e)
         enddo
       enddo
     enddo
!$OMP END PARALLEL DO

     write(6,'(i2,x,a10,x,16(l,x))') crystallite_Niteration,'cryst_onT',crystallite_onTrack
 write (6,'(a,/,3(3(f12.7,x)/))') 'Lp of 1 8 1',crystallite_Lp(1:3,:,1,8,1)
     
! --+>> state integration <<+--
!
! incrementing by crystallite_subdt
! based on constitutive_subState0
! results in constitutive_state

!$OMP PARALLEL DO
     do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
         do g = 1,myNgrains
           if (crystallite_requested(g,i,e) .and. &
               crystallite_onTrack(g,i,e)) then                   ! all undone crystallites
             crystallite_converged(g,i,e) = crystallite_updateState(g,i,e)
             if (crystallite_converged(g,i,e)) then
!$OMP CRITICAL (distributionState)
               debug_StateLoopDistribution(crystallite_Niteration) = &
               debug_StateLoopDistribution(crystallite_Niteration) + 1
!$OMP END CRITICAL (distributionState)
             endif
           endif
         enddo
       enddo
     enddo
!$OMP END PARALLEL DO

   enddo                                                         ! crystallite convergence loop  

   write(6,*)
   write(6,'(a10,x,16(f6.4,x))') 'cryst_frac',crystallite_subFrac
   write(6,'(a10,x,16(f6.4,x))') 'cryst_step',crystallite_subStep
   write(6,'(a10,x,16(l,x))') 'cryst_req',crystallite_requested
   write(6,'(a10,x,16(l,x))') 'cryst_onT',crystallite_onTrack
   write(6,'(a10,x,16(l,x))') 'cryst_cvg',crystallite_converged
   write(6,'(a10,x,16(e8.3,x))') 'cryst_dt',crystallite_subdt
 enddo                                                           ! cutback loop

 

! ------ check for non-converged crystallites ------

!$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
     do g = 1,myNgrains
       if (.not. crystallite_converged(g,i,e)) then           ! respond fully elastically
         call IO_warning(600,e,i,g)
         invFp = math_inv3x3(crystallite_partionedFp0(:,:,g,i,e))
         Fe_guess = math_mul33x33(crystallite_partionedF(:,:,g,i,e),invFp)
         PK2 = math_Mandel6to33( &
                 math_mul66x6( &
                   0.5_pReal*constitutive_homogenizedC(g,i,e), &
                   math_Mandel33to6( &
                     math_mul33x33(transpose(Fe_guess),Fe_guess) - math_I3 &
                   ) &
                 ) &
               )
         crystallite_P(:,:,g,i,e) = math_mul33x33(Fe_guess,math_mul33x33(PK2,transpose(invFp)))
       endif
     enddo
   enddo
 enddo
!$OMP END PARALLEL DO

! ------ stiffness calculation ------

 if(updateJaco) then                                                  ! Jacobian required
   if (debugger) then
!$OMP CRITICAL (write2out)
     write (6,*) 'Jacobian calc'
     write(6,'(a10,x,16(f6.4,x))') 'cryst_dt',crystallite_subdt
!$OMP END CRITICAL (write2out)
   endif

!$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
       do g = 1,myNgrains
         if (crystallite_converged(g,i,e)) then                         ! grain converged in above iteration
           mySizeState = constitutive_sizeState(g,i,e)                  ! number of state variables for this grain
           myState(1:mySizeState) = constitutive_state(g,i,e)%p         ! remember unperturbed, converged state...
           myF  = crystallite_subF(:,:,g,i,e)                           ! ... and kinematics
           myFp = crystallite_Fp(:,:,g,i,e)
           myFe = crystallite_Fe(:,:,g,i,e)
           myLp = crystallite_Lp(:,:,g,i,e)
           myP  = crystallite_P(:,:,g,i,e)
           do k = 1,3                                                   ! perturbation...
              do l = 1,3                                                 ! ...components
                crystallite_subF(:,:,g,i,e) =  myF                        ! initialize perturbed F to match converged
                crystallite_subF(k,l,g,i,e) = crystallite_subF(k,l,g,i,e) + pert_Fg  ! perturb single component
                converged = .false.
                crystallite_Niteration = 0_pInt
                do while(.not. converged .and. crystallite_Niteration < nCryst) ! keep cycling until done (potentially non-converged)
                  crystallite_Niteration = crystallite_Niteration + 1_pInt
                  if(crystallite_integrateStress(g,i,e))  &               ! stress of perturbed situation (overwrites_P,_Tstar_v,_Fp,_Lp,_Fe)
                    converged = crystallite_updateState(g,i,e)
                end do
                if (converged)  &                                         ! converged state warrants  stiffness update
                  crystallite_dPdF(:,:,k,l,g,i,e) =(crystallite_p(:,:,g,i,e) - myP)/pert_Fg ! tangent dP_ij/dFg_kl 
!$OMP CRITICAL (out)
                debug_StiffnessStateLoopDistribution(crystallite_Niteration) = &
                  debug_StiffnessstateLoopDistribution(crystallite_Niteration) + 1 
!$OMP END CRITICAL (out)
              end do
            end do
           constitutive_state(g,i,e)%p = myState                        ! restore unperturbed, converged state...
           crystallite_Fp(:,:,g,i,e) = myFp                             ! ... and kinematics
           crystallite_Fe(:,:,g,i,e) = myFe
           crystallite_Lp(:,:,g,i,e) = myLp
           crystallite_P(:,:,g,i,e)  = myP
         else                                                           ! grain has not converged
           crystallite_dPdF(:,:,:,:,g,i,e) = crystallite_fallbackdPdF(:,:,:,:,g,i,e)  ! use fallback
         endif
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO
 endif
 
end subroutine




!********************************************************************
! update the internal state of the constitutive law
! and tell whether state has converged
!********************************************************************
 function crystallite_updateState(&
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )
 use prec, only: pReal,pInt,rTol_crystalliteState
 use constitutive, only: constitutive_dotState,constitutive_sizeDotState,&
                         constitutive_subState0,constitutive_state
 use debug
 
 logical crystallite_updateState

 integer(pLongInt) tick,tock,tickrate,maxticks
 integer(pInt) g,i,e,mySize
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(constitutive_sizeDotState(g,i,e)) :: residuum

 mySize = constitutive_sizeDotState(g,i,e)
 call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 residuum = constitutive_state(g,i,e)%p(1:mySize) - constitutive_subState0(g,i,e)%p(1:mySize) - &
            crystallite_subdt(g,i,e)*&
            constitutive_dotState(crystallite_Tstar_v(:,g,i,e),crystallite_Temperature(g,i,e),g,i,e)   ! residuum from evolution of microstructure
 call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
 debug_cumDotStateCalls = debug_cumDotStateCalls + 1_pInt
 debug_cumDotStateTicks = debug_cumDotStateTicks + tock-tick
 if (tock < tick) debug_cumDotStateTicks  = debug_cumDotStateTicks + maxticks
 constitutive_state(g,i,e)%p(1:mySize) = constitutive_state(g,i,e)%p(1:mySize) - residuum              ! update of microstructure
 crystallite_updateState = maxval(abs(residuum/constitutive_state(g,i,e)%p(1:mySize)),&
                                  constitutive_state(g,i,e)%p(1:mySize) /= 0.0_pReal) < rTol_crystalliteState
 return

 end function



!***********************************************************************
!***     calculation of stress (P), stiffness (dPdF),                ***
!***     and announcement of any                                     ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 function crystallite_integrateStress(&
     g,&          ! grain number
     i,&          ! integration point number
     e)           ! element number

 use prec, only: pReal,pInt,subStepMin
 use debug
 use constitutive, only: constitutive_state,constitutive_subState0,constitutive_sizeState
 use math

 implicit none

 integer(pInt), intent(in) :: e,i,g
 integer(pInt) mySize, Niteration
 real(pReal) dt_aim,subFrac,subStep,det
 real(pReal), dimension(3,3)     :: inv
 real(pReal), dimension(3,3)     :: Fg_current,Fg_aim,deltaFg
 real(pReal), dimension(3,3)     :: Fp_current,Fp_new
 real(pReal), dimension(3,3)     :: Fe_current,Fe_new
 real(pReal), dimension(constitutive_sizeState(g,i,e)) :: interpolatedState
 
 logical crystallite_integrateStress              ! still on track?

 mySize = constitutive_sizeState(g,i,e)
 deltaFg = crystallite_subF(:,:,g,i,e) - crystallite_subF0(:,:,g,i,e)
 Fg_current = crystallite_subF0(:,:,g,i,e)                        ! initialize to start of inc
 Fp_current = crystallite_subFp0(:,:,g,i,e)
 Fe_current = math_mul33x33(Fg_current,math_inv3x3(Fp_current))
 subFrac = 0.0_pReal
 subStep = 1.0_pReal
 Niteration = 0_pInt
 crystallite_integrateStress = .false.                             ! be pessimisitc

! begin the cutback loop
 do while (subStep > subStepMin)                                   ! continue until finished or too much cut backing
   Niteration = Niteration + 1
   Fg_aim = Fg_current + subStep*deltaFg                           ! aim for Fg
   dt_aim = subStep*crystallite_subdt(g,i,e)                       ! aim for dt
   debugger = (g == 1 .and. i == 8 .and. e == 1)
   call TimeIntegration(crystallite_integrateStress,&
                        crystallite_Lp(:,:,g,i,e),crystallite_Fp(:,:,g,i,e),crystallite_Fe(:,:,g,i,e),&
                        crystallite_Tstar_v(:,g,i,e),crystallite_P(:,:,g,i,e), &
                        Fg_aim,Fp_current,crystallite_Temperature(g,i,e),&
                        (          subFrac+subStep)*constitutive_state(g,i,e)%p(1:mySize) + &
                        (1.0_pReal-subFrac-subStep)*constitutive_subState0(g,i,e)%p(1:mySize),&    ! interpolated state
                        dt_aim,g,i,e)
   if (crystallite_integrateStress) then                           ! happy with time integration
     if (e == 1 .and. i == 8 .and. g == 1) then
       write(6,*) '*** winding forward in IntegrateStress ***'
       write(6,*) subFrac, subFrac+subStep
       write (6,'(a,/,3(3(f12.7,x)/))') 'Lp of 1 8 1',crystallite_Lp(1:3,:,1,8,1)
       write (6,'(a,/,3(3(f12.7,x)/))') 'Fp of 1 8 1',crystallite_Fp(1:3,:,1,8,1)
     endif
     Fg_current = Fg_aim                                           ! wind forward
     Fe_current = crystallite_Fe(:,:,g,i,e)
     Fp_current = crystallite_Fp(:,:,g,i,e)
     subFrac = subFrac + subStep
     subStep = min(1.0_pReal-subFrac, subStep*2.0_pReal)           ! accelerate
   else                                                            ! time integration encountered trouble
     subStep = 0.5_pReal * subStep                                 ! cut time step in half
     crystallite_Lp(:,:,g,i,e) = 0.5_pReal*(crystallite_Lp(:,:,g,i,e) + &  ! interpolate Lp  and L
                                            (math_I3 - math_mul33x33(Fp_current,&
                                                       math_mul33x33(math_inv3x3(Fg_aim),Fe_current)))/dt_aim)
   endif
 enddo  ! potential substepping

!$OMP CRITICAL (distributionStress)
 debug_StressLoopDistribution(Niteration) = debug_StressLoopDistribution(Niteration) + 1
!$OMP END CRITICAL (distributionStress)

 return                                                            ! with final happyness

 end function


!***********************************************************************
!***     fully-implicit two-level time integration                   ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
subroutine TimeIntegration(&
     happy,&            ! return status
     Lpguess,&          ! guess of plastic velocity gradient
     Fp_new,&           ! new plastic deformation gradient
     Fe_new,&           ! new "elastic" deformation gradient
     Tstar_v,&          ! Stress vector
     P,&                ! 1st PK stress (taken as initial guess if /= 0)
!
     Fg_new,&           ! new total def gradient
     Fp_old,&           ! former plastic def gradient
     Temperature,&      ! temperature
     state,&            ! microstructural state
     dt,&               ! time increment
     grain,&            ! grain number
     ip,&               ! integration point number
     cp_en &            ! element number
    )
    
 use prec
 use debug
 use mesh, only: mesh_element
 use constitutive, only: constitutive_microstructure,constitutive_homogenizedC,constitutive_LpAndItsTangent,&
                         constitutive_sizeState
 use math
 use IO
 implicit none

 logical, intent(out) :: happy
 real(pReal), dimension(3,3), intent(inout) :: Lpguess
 real(pReal), dimension(3,3), intent(out) :: Fp_new,Fe_new,P
 real(pReal), dimension(6),   intent(out) :: Tstar_v
 real(pReal), dimension(3,3), intent(in)  :: Fg_new,Fp_old
 real(pReal), intent(in) :: Temperature,dt
 integer(pInt), intent(in) :: cp_en, ip, grain
 real(pReal), dimension(constitutive_sizeState(grain,ip,cp_en)), intent(in) :: state

 logical error
 integer(pInt) Niteration,dummy, i,j,k,l,m,n
 integer(pLongInt) tick,tock,tickrate,maxticks
 real(pReal) p_hydro,det, leapfrog,maxleap
 real(pReal), dimension(3,3,3,3) :: C
 real(pReal), dimension(9,9) :: dLp,dTdLp,dRdLp,invdRdLp,eye2
 real(pReal), dimension(6,6) :: C_66
 real(pReal), dimension(3,3) :: invFp_new,invFp_old,Lp,Lpguess_old,Rinner,Rinner_old,A,B,BT,AB,BTA

 happy = .false.                            ! be pessimistic
 eye2 = math_identity2nd(9)

 invFp_old = math_inv3x3(Fp_old)            ! inversion of Fp_old
 if (all(invFp_old == 0.0_pReal)) return    ! failed
! write (6,'(a,/,3(3(f12.7,x)/))') 'Fp old',Fp_old
! write (6,'(a,/,3(3(f12.7,x)/))') 'Fp old inv',invFp_old
 
 A = math_mul33x33(transpose(invFp_old), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_old)))

!$OMP CRITICAL (write2out)
! if (debugger) write (6,'(a,/,3(3(f12.7,x)/))') 'Fg to be calculated',Fg_new
!$OMP END CRITICAL (write2out)

 call constitutive_microstructure(Temperature,grain,ip,cp_en)
 C_66 = constitutive_homogenizedC(grain,ip,cp_en)
 C = math_Mandel66to3333(C_66)         ! 4th rank elasticity tensor

 Niteration = 0_pInt
 leapfrog = 1.0_pReal                  ! correction as suggested by invdRdLp-step
 maxleap = 1024.0_pReal                ! preassign maximum acceleration level

 Lpguess_old = Lpguess                 ! consider present Lpguess good (i.e. worth remembering)

Inner: do              ! inner iteration: Lp
   Niteration = Niteration + 1
   if (Niteration > nLp) then          ! too many loops required
     Lpguess = Lpguess_old             ! do not trust the last update but resort to former one
     return
   endif

!   write(6,*) 'iteration',Niteration
!   write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess

   B = math_i3 - dt*Lpguess
   BT = transpose(B)
   AB = math_mul33x33(A,B)
   BTA = math_mul33x33(BT,A)
   Tstar_v = 0.5_pReal*math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB)-math_I3))
   p_hydro = (Tstar_v(1) + Tstar_v(2) + Tstar_v(3))/3.0_pReal
   forall(i=1:3) Tstar_v(i) = Tstar_v(i) - p_hydro              ! subtract hydrostatic pressure
!   write (6,'(a,/,6(f12.7,x))') 'Tstar',Tstar_v
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
   call constitutive_LpAndItsTangent(Lp,dLp, Tstar_v,Temperature,grain,ip,cp_en)
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   debug_cumLpCalls = debug_cumLpCalls + 1_pInt
   debug_cumLpTicks  = debug_cumLpTicks + tock-tick
   if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks

   Rinner = Lpguess - Lp                                        ! update current residuum
!   write (6,'(a,/,3(3(f12.7,x)/))') 'Lp',Lp
!   write (6,'(a,/,3(3(f12.7,x)/))') 'Residuum',Rinner

   if (.not.(any(Rinner/=Rinner)) .and. &                       ! exclude any NaN in residuum
       ( ( maxval(abs(Rinner)) < aTol_crystalliteStress) .or. & ! below abs tol .or.
         ( any(abs(dt*Lpguess) > relevantStrain) .and. &        ! worth checking? .and.
           maxval(abs(Rinner/Lpguess),abs(dt*Lpguess) > relevantStrain) < rTol_crystalliteStress &  ! below rel tol
         ) &
       )   &
      )    &
     exit Inner                                                 ! convergence
!
!          check for acceleration/deceleration in Newton--Raphson correction
!
   if (any(Rinner/=Rinner) .and. &                              ! NaN occured at regular speed
       leapfrog == 1.0) then
     Lpguess = Lpguess_old                                      ! restore known good guess and croak for cutback
     return

   elseif (leapfrog > 1.0_pReal .and. &                         ! at fast pace ?
           (sum(Rinner*Rinner) > sum(Rinner_old*Rinner_old) .or. &  ! worse residuum
            sum(Rinner*Rinner_old) < 0.0_pReal) .or. &              ! residuum changed sign (overshoot)
           any(Rinner/=Rinner) ) then                           ! NaN
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
     call math_invert(9,dRdLp,invdRdLp,dummy,error)             ! invert dR/dLp --> dLp/dR
     if (error) then
       if (debugger) then
!$OMP CRITICAL (write2out)
                 write (6,*) 'inversion dR/dLp failed',grain,ip,cp_en
!                 write (6,'(a,/,9(9(e9.3,x)/))') 'dRdLp', dRdLp(1:9,:)
!                 write (6,'(a,/,3(4(f9.3,x)/))') 'state_new / MPa',state/1e6_pReal
                 write (6,'(a,/,3(3(f12.7,x)/))') 'Lpguess',Lpguess(1:3,:)
                 write (6,'(a,/,3(3(e12.7,x)/))') 'Lp',Lp(1:3,:)
                 write (6,'(a,/,6(f9.3,x))') 'Tstar / MPa',Tstar_v/1e6_pReal
!$OMP END CRITICAL (write2out)
       endif 
       return
     endif

     Rinner_old = Rinner                                        ! remember current residuum
     Lpguess_old = Lpguess                                      ! remember current Lp guess 
     if (Niteration > 1 .and. leapfrog < maxleap) leapfrog = 2.0_pReal * leapfrog   ! accelerate if ok
   endif

   Lpguess = Lpguess_old                                        ! start from current guess                                       
   Rinner  = Rinner_old                                         ! use current residuum
   forall (i=1:3,j=1:3,k=1:3,l=1:3) &                           ! leapfrog to updated Lpguess 
     Lpguess(i,j) = Lpguess(i,j) - leapfrog*invdRdLp(3*(i-1)+j,3*(k-1)+l)*Rinner(k,l)
 enddo Inner

!$OMP CRITICAL (distributionLp)
 debug_LpLoopDistribution(Niteration) = debug_LpLoopDistribution(Niteration) + 1
!$OMP END CRITICAL (distributionLp)
 invFp_new = math_mul33x33(invFp_old,B)
 if (debugger) then
   write (6,'(a,x,f10.6,/,3(3(f12.7,x)/))') 'Lp(guess)',dt,Lpguess(1:3,:)
   write (6,'(a,/,3(3(f12.7,x)/))') 'invFp_old',invFp_old(1:3,:)
   write (6,'(a,/,3(3(f12.7,x)/))') 'B',B(1:3,:)
   write (6,'(a,/,3(3(f12.7,x)/))') 'invFp_new',invFp_new(1:3,:)
 endif
 call math_invert3x3(invFp_new,Fp_new,det,error)
 if (debugger) then
   write (6,'(a,/,3(3(f12.7,x)/))') 'invFp_new',invFp_new(1:3,:)
   write (6,'(a,/,3(3(f12.7,x)/))') 'Fp_new',Fp_new(1:3,:)
   write (6,'(a,x,l,x,a,f10.6)') 'with inversion error:',error,'and determinant:',det
 endif
 if (error) return                                             ! inversion failed

 Fp_new = Fp_new*det**(1.0_pReal/3.0_pReal)                    ! regularize Fp by det = det(InvFp_new) !!
 forall (i=1:3) Tstar_v(i) = Tstar_v(i) + p_hydro              ! add hydrostatic component back
 Fe_new = math_mul33x33(Fg_new,invFp_new)                      ! calc resulting Fe
 P = math_mul33x33(Fe_new,math_mul33x33(math_Mandel6to33(Tstar_v),transpose(invFp_new)))    ! first PK stress

 happy = .true.                                                ! now smile...

 return

end subroutine


 
!********************************************************************
! return results of particular grain
!********************************************************************
function crystallite_postResults(&
   Tstar_v,&        ! stress
   Temperature, &   ! temperature
   dt,&             ! time increment
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )

 use prec,     only: pInt,pReal
 use math,     only: math_pDecomposition,math_RtoEuler, inDeg
 use IO,       only: IO_warning
 use material, only: material_phase,material_volume
 use constitutive, only: constitutive_sizePostResults, constitutive_postResults
 implicit none

 integer(pInt), intent(in) :: g,i,e
 real(pReal), intent(in) ::   Temperature,dt
 real(pReal), dimension(6), intent(in) :: Tstar_v
 real(pReal), dimension(3,3) :: U,R
 logical error

 real(pReal), dimension(crystallite_Nresults + constitutive_sizePostResults(g,i,e)) :: crystallite_postResults
 
 if (crystallite_Nresults >= 2) then
   crystallite_postResults(1) = material_phase(g,i,e)
   crystallite_postResults(2) = material_volume(g,i,e)
 endif
 if (crystallite_Nresults >= 5) then
   call math_pDecomposition(crystallite_Fe(:,:,g,i,e),U,R,error)               ! polar decomposition of Fe
   if (error) then
     call IO_warning(650,e,i,g)
     crystallite_postResults(3:5) = (/400.0,400.0,400.0/)                 ! fake orientation
   else
     crystallite_postResults(3:5) = math_RtoEuler(transpose(R))*inDeg     ! orientation
   endif
 endif
 
 crystallite_postResults(crystallite_Nresults+1:crystallite_Nresults+constitutive_sizePostResults(g,i,e)) = &
   constitutive_postResults(Tstar_v,Temperature,dt,g,i,e)
 return 
 
end function


END MODULE
!##############################################################