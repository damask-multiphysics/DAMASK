
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

 endsubroutine


!********************************************************************
! calculate stress (P) and tangent (dPdF) for crystallites
!********************************************************************
subroutine crystallite_stressAndItsTangent(updateJaco)

 use prec,         only: pInt,pReal,subStepMin,nCryst
 use debug
 use IO,           only: IO_warning
 use math
 use FEsolving,    only: FEsolving_execElem, FEsolving_execIP, theInc
 use mesh,         only: mesh_element
 use material,     only: homogenization_Ngrains
 use constitutive
 implicit none
 
 logical, intent(in) :: updateJaco
 real(pReal), dimension(3,3) :: invFp,Fe_guess,PK2,myF,myFp,myFe,myLp,myP
 real(pReal), dimension(constitutive_maxSizeState) :: myState
 integer(pInt) NiterationCrystallite, NiterationState
 integer(pInt) g,i,e,k,l, myNgrains, mySizeState
 logical onTrack,converged

! ------ initialize to starting condition ------

 write (6,*)
 write (6,*) 'Crystallite request from Materialpoint'
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedF0  of 1 1 1',crystallite_partionedF0(1:3,:,1,1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedFp0 of 1 1 1',crystallite_partionedFp0(1:3,:,1,1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedF   of 1 1 1',crystallite_partionedF(1:3,:,1,1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'crystallite_partionedLp0 of 1 1 1',crystallite_partionedLp0(1:3,:,1,1,1)


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

 NiterationCrystallite = 0_pInt

 do while (any(crystallite_subStep(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMin))  ! cutback loop for crystallites
   
   NiterationCrystallite = NiterationCrystallite + 1
   
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
         debugger = (g == 1 .and. i == 1 .and. e == 1)
         if (crystallite_converged(g,i,e)) then
           crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
           crystallite_subStep(g,i,e) = min(1.0_pReal-crystallite_subFrac(g,i,e), 2.0_pReal * crystallite_subStep(g,i,e))
           if (crystallite_subStep(g,i,e) > subStepMin) then
             crystallite_subF0(:,:,g,i,e)    = crystallite_subF(:,:,g,i,e)         ! wind forward...
             crystallite_subFp0(:,:,g,i,e)   = crystallite_Fp(:,:,g,i,e)           ! ...plastic def grad
             crystallite_subLp0(:,:,g,i,e)   = crystallite_Lp(:,:,g,i,e)           ! ...plastic velocity gradient
             constitutive_subState0(g,i,e)%p = constitutive_state(g,i,e)%p         ! ...microstructure
           endif
           if (debugger) then
!$OMP CRITICAL (write2out)
             write(6,'(a21,f6.4,a28,f6.4,a35)') 'winding forward from ', &
                crystallite_subFrac(g,i,e)-crystallite_subStep(g,i,e),' to new crystallite_subfrac ', &
                crystallite_subFrac(g,i,e),' in crystallite_stressAndItsTangent'
             write(6,*)
!$OMP END CRITICAL (write2out)
           endif
         else
           crystallite_subStep(g,i,e) = 0.5_pReal * crystallite_subStep(g,i,e)     ! cut step in half and restore...
           crystallite_Fp(:,:,g,i,e) = crystallite_subFp0(:,:,g,i,e)               ! ...plastic def grad
           crystallite_Lp(:,:,g,i,e) = crystallite_subLp0(:,:,g,i,e)               ! ...plastic velocity grad
           constitutive_state(g,i,e)%p = constitutive_subState0(g,i,e)%p           ! ...microstructure
           if (debugger) then
!$OMP CRITICAL (write2out)
             write(6,'(a78,f6.4)') 'cutback step in crystallite_stressAndItsTangent with new crystallite_subStep: ',&
                                    crystallite_subStep(g,i,e)
             write(6,*)
!$OMP END CRITICAL (write2out)
           endif
         endif
     
         crystallite_onTrack(g,i,e) = crystallite_subStep(g,i,e) > subStepMin      ! still on track or already done (beyond repair)
         if (crystallite_onTrack(g,i,e)) then                                      ! specify task (according to substep)
           crystallite_subF(:,:,g,i,e)  = crystallite_subF0(:,:,g,i,e) + &
                                          crystallite_subStep(g,i,e) * &
                                          (crystallite_partionedF(:,:,g,i,e) - crystallite_partionedF0(:,:,g,i,e))
           crystallite_subdt(g,i,e)     = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
           if (debugger) then
!$OMP CRITICAL (write2out)
             write(6,'(a36,e8.3)') 'current timestep crystallite_subdt: ',crystallite_subdt(g,i,e)
             write(6,*)
!$OMP END CRITICAL (write2out)
           endif
           crystallite_converged(g,i,e) = .false.                                  ! start out non-converged
         endif
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO

! ------ convergence loop for stress and state ------

   NiterationState = 0_pInt
   if (debugger)  write(6,*) 'state integration started'
   
   do while (any(            crystallite_requested(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and.       crystallite_onTrack(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. crystallite_converged(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. NiterationState < nState)              ! convergence loop for crystallite
     
     NiterationState = NiterationState + 1

! --+>> stress integration <<+--
!
! incrementing by crystallite_subdt
! based on crystallite_subF0,.._subFp0,.._subLp0
!          constitutive_state is internally interpolated with .._subState0
!          to account for substepping within _integrateStress
! results in crystallite_Fp,.._Lp

     if (debugger)  write(6,*) 'stress integration started'

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

     if (crystallite_requested(1,1,1) .and. crystallite_onTrack(1,1,1)) then
       write(6,*) 'stress integration converged'
       write(6,'(a,/,3(3(e15.7,x)/))') 'P of 1 1 1',crystallite_P(:,:,1,1,1)
       write(6,'(a,/,3(3(f12.7,x)/))') 'Lp of 1 1 1',crystallite_Lp(:,:,1,1,1)
     endif
     
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
               debug_StateLoopDistribution(NiterationState) = debug_StateLoopDistribution(NiterationState) + 1
!$OMP END CRITICAL (distributionState)
!$OMP CRITICAL (distributionCrystallite)
               debug_CrystalliteLoopDistribution(NiterationCrystallite) = &
                  debug_CrystalliteLoopDistribution(NiterationCrystallite) + 1
!$OMP END CRITICAL (distributionCrystallite)
             endif
           endif
         enddo
       enddo
     enddo
!$OMP END PARALLEL DO

     if (crystallite_requested(1,1,1) .and. crystallite_onTrack(1,1,1) .and. crystallite_converged(1,1,1)) then
       write(6,*) 'state integration converged'
       write(6,'(a20,e8.3)') 'state of 1 1 1: ', constitutive_state(1,1,1)%p(1)
       write(6,*)
     endif
     
   enddo                                                         ! crystallite convergence loop  

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
     write (6,*) 'Stiffness calculation started'
!$OMP CRITICAL (write2out)  
!     write(6,'(a10,x,16(f6.4,x))') 'cryst_dt',crystallite_subdt
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
		   if (debugger) then
			 write (6,*) '#############'
			 write (6,*) 'central solution'
			 write (6,*) '#############'
			 write (6,'(a,/,3(3(f12.4,x)/))') '    P of 1 1 1',myP(1:3,:)/1e6
			 write (6,'(a,/,3(3(f12.8,x)/))') '   Fp of 1 1 1',myFp(1:3,:)
			 write (6,'(a,/,3(3(f12.8,x)/))') '   Lp of 1 1 1',myLp(1:3,:)
			 write (6,'(a,/,f12.4)') 'state of 1 1 1',myState/1e6
		   endif
           do k = 1,3                                                   ! perturbation...
             do l = 1,3                                                 ! ...components
               crystallite_subF(:,:,g,i,e) = myF                        ! initialize perturbed F to match converged
               crystallite_subF(k,l,g,i,e) = crystallite_subF(k,l,g,i,e) + pert_Fg  ! perturb single component
               onTrack = .true.
               converged = .false.
               NiterationState = 0_pInt
               if (debugger) then
                 write (6,*) '============='
                 write (6,'(i1,x,i1)') k,l
                 write (6,*) '============='
                 write (6,'(a,/,3(3(f12.6,x)/))') 'pertF of 1 1 1',crystallite_subF(1:3,:,g,i,e)
               endif
               do while(.not. converged .and. onTrack .and. NiterationState < nState) ! keep cycling until done (potentially non-converged)
                 NiterationState = NiterationState + 1_pInt
                 if (debugger) write (6,'(a4,x,i6)') 'loop',NiterationState
                 onTrack = crystallite_integrateStress(g,i,e)           ! stress of perturbed situation (overwrites _P,_Tstar_v,_Fp,_Lp,_Fe)
                 if(onTrack) converged = crystallite_updateState(g,i,e)
                 if (debugger) then
                 write (6,*) '-------------'
                   write (6,'(l,x,l)') onTrack,converged
                   write (6,'(a,/,3(3(f12.4,x)/))') 'pertP of 1 1 1',crystallite_P(1:3,:,g,i,e)/1e6
                   write (6,'(a,/,3(3(f12.4,x)/))') 'DP    of 1 1 1',(crystallite_P(1:3,:,g,i,e)-myP(1:3,:))/1e6
                   write (6,'(a,/,f12.4)') 'state  of 1 1 1',constitutive_state(g,i,e)%p/1e6
                   write (6,'(a,/,f12.4)') 'Dstate of 1 1 1',(constitutive_state(g,i,e)%p-myState)/1e6
                 endif
               enddo
               if (converged) &                                         ! converged state warrants stiffness update
                 crystallite_dPdF(:,:,k,l,g,i,e) = (crystallite_P(:,:,g,i,e) - myP)/pert_Fg ! tangent dP_ij/dFg_kl
               constitutive_state(g,i,e)%p = myState                        ! restore unperturbed, converged state...
               crystallite_Fp(:,:,g,i,e) = myFp                             ! ... and kinematics
               crystallite_Fe(:,:,g,i,e) = myFe
               crystallite_Lp(:,:,g,i,e) = myLp
               crystallite_P(:,:,g,i,e)  = myP
!$OMP CRITICAL (out)
               debug_StiffnessStateLoopDistribution(NiterationState) = &
               debug_StiffnessstateLoopDistribution(NiterationState) + 1
!$OMP END CRITICAL (out)
             enddo
           enddo
           constitutive_state(g,i,e)%p = myState                        ! restore unperturbed, converged state...
           crystallite_Fp(:,:,g,i,e) = myFp                             ! ... and kinematics
           crystallite_Fe(:,:,g,i,e) = myFe
           crystallite_Lp(:,:,g,i,e) = myLp
           crystallite_P(:,:,g,i,e)  = myP
           if (e == 1 .and. i == 1 .and. g == 1) then
             write (6,'(a,/,9(9(f12.4,x)/))') 'dPdF/GPa',crystallite_dPdF(:,:,:,:,1,1,1)/1e9
           endif
         else                                                           ! grain has not converged
           crystallite_dPdF(:,:,:,:,g,i,e) = crystallite_fallbackdPdF(:,:,:,:,g,i,e)  ! use fallback
         endif
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO
 endif
 
endsubroutine




!********************************************************************
! update the internal state of the constitutive law
! and tell whether state has converged
!********************************************************************
 function crystallite_updateState(&
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )
 use prec, only:  pReal, &
                  pInt, &
                  rTol_crystalliteState
 use constitutive, only:  constitutive_dotState, &
                          constitutive_sizeDotState, &
                          constitutive_subState0, &
                          constitutive_state
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
 
 if (any(constitutive_state(g,i,e)%p(1:mySize)/=constitutive_state(g,i,e)%p(1:mySize))) return           ! NaN occured?
 
 constitutive_state(g,i,e)%p(1:mySize) = constitutive_state(g,i,e)%p(1:mySize) - residuum              ! update of microstructure
 crystallite_updateState = maxval(abs(residuum/constitutive_state(g,i,e)%p(1:mySize)),&
                                  constitutive_state(g,i,e)%p(1:mySize) /= 0.0_pReal) < rTol_crystalliteState
 return

 endfunction



!***********************************************************************
!***     calculation of stress (P) with time integration             ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 function crystallite_integrateStress(&
     g,&          ! grain number
     i,&          ! integration point number
     e)           ! element number

 !*** variables and functions from other modules ***!
 use prec, only:                      pReal, &
                                      pInt, &
                                      nStress, &
                                      aTol_crystalliteStress, &
                                      rTol_crystalliteStress, &
                                      relevantStrain
 use debug, only:                     debugger, &
                                      debug_cumLpCalls, &
                                      debug_cumLpTicks, &
                                      debug_StressLoopDistribution
 use constitutive, only:              constitutive_microstructure, &
                                      constitutive_homogenizedC, &
                                      constitutive_LpAndItsTangent
 use math, only:                      math_mul33x33, &
                                      math_mul66x6, &
                                      math_mul99x99, &
                                      math_inv3x3, &
                                      math_invert3x3, &
                                      math_i3, &
                                      math_identity2nd, &
                                      math_Mandel66to3333, &
                                      math_mandel33to6

 implicit none

 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 
 !*** output variables ***!
 logical                              crystallite_integrateStress   ! flag indicating if integration suceeded
 
 !*** internal local variables ***!
 real(pReal), dimension(3,3)::        Fg_current, &                 ! deformation gradient at start of timestep
                                      Fg_new, &                     ! deformation gradient at end of timestep
                                      Fp_current, &                 ! plastic deformation gradient at start of timestep
                                      Fp_new, &                     ! plastic deformation gradient at end of timestep
                                      Fe_current, &                 ! elastic deformation gradient at start of timestep
                                      Fe_new, &                     ! elastic deformation gradient at end of timestep
                                      invFp_new, &                  ! inverse of Fp_new
                                      invFp_current, &              ! inverse of Fp_current
                                      Lpguess, &                    ! current guess for plastic velocity gradient
                                      Lpguess_old, &                ! known last good guess for plastic velocity gradient
                                      Lp_constitutive, &            ! plastic velocity gradient resulting from constitutive law
                                      residuum, &                   ! current residuum of plastic velocity gradient
                                      residuum_old, &               ! last residuum of plastic velocity gradient
                                      A, &
                                      B, &
                                      BT, &
                                      AB, &
                                      BTA
 real(pReal), dimension(6)::          Tstar_v                       ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
 real(pReal), dimension(9,9)::        dLp_constitutive, &           ! partial derivative of plastic velocity gradient calculated by constitutive law
                                      dTdLp, &                      ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dRdLp, &                      ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                      invdRdLp                      ! inverse of dRdLp
 real(pReal), dimension(3,3,3,3)::    C                             ! 4th rank elasticity tensor
 real(pReal), dimension(6,6)::        C_66                          ! simplified 2nd rank elasticity tensor 
 real(pReal)                          p_hydro, &                    ! volumetric part of 2nd Piola-Kirchhoff Stress
                                      det, &                        ! determinant
                                      leapfrog, &                   ! acceleration factor for Newton-Raphson scheme
                                      maxleap                       ! maximum acceleration factor
 logical                              error                         ! flag indicating an error
 integer(pInt)                        NiterationStress, &           ! number of stress integrations
                                      dummy, &
                                      h, &
                                      j, &
                                      k, &
                                      l, &
                                      m, &
                                      n
 integer(pLongInt)                    tick, &
                                      tock, &
                                      tickrate, &
                                      maxticks
 
 !*** global variables ***!
 ! crystallite_subF0
 ! crystallite_subF
 ! crystallite_subFp0
 ! crystallite_Tstar_v
 ! crystallite_Lp
 ! crystallite_subdt
 ! crystallite_Temperature
 

 ! be pessimistic
 crystallite_integrateStress = .false.

 ! feed local variables
 Fg_current =   crystallite_subF0(:,:,g,i,e)
 Fg_new =       crystallite_subF(:,:,g,i,e)
 Fp_current =   crystallite_subFp0(:,:,g,i,e)
 Fe_current =   math_mul33x33(Fg_current,math_inv3x3(Fp_current))
 Tstar_v =      crystallite_Tstar_v(:,g,i,e)
 Lpguess_old =  crystallite_Lp(:,:,g,i,e)                           ! consider present Lp good (i.e. worth remembering) ...
 Lpguess =      crystallite_Lp(:,:,g,i,e)                           ! ... and take it as first guess
 
 ! inversion of Fp_current...
 invFp_current = math_inv3x3(Fp_current)                            
 if (all(invFp_current == 0.0_pReal)) then                          ! ... failed?
   if (debugger) write(6,*) '::: integrateStress failed on invFp_current inversion'
   return
 endif
 
 A = math_mul33x33(transpose(invFp_current), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_current)))
 
 ! update microstructure
 call constitutive_microstructure(crystallite_Temperature(g,i,e),g,i,e)
 
 ! get elasticity tensor
 C_66 = constitutive_homogenizedC(g,i,e)
 C = math_Mandel66to3333(C_66)
 
 ! start LpLoop with no acceleration
 NiterationStress = 0_pInt
 leapfrog = 1.0_pReal
 maxleap = 1024.0_pReal

LpLoop: do
   
   ! increase loop counter
   NiterationStress = NiterationStress + 1
   
   ! too many loops required ?
   if (NiterationStress > nStress) then
     if (debugger) write(6,*) '::: integrateStress exceeded nStress loopcount'
     return
   endif
   
   B = math_i3 - crystallite_subdt(g,i,e)*Lpguess
   BT = transpose(B)
   AB = math_mul33x33(A,B)
   BTA = math_mul33x33(BT,A)
   
   ! calculate 2nd Piola-Kirchhoff stress tensor
   Tstar_v = 0.5_pReal*math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB)-math_I3))
   p_hydro = sum(Tstar_v(1:3))/3.0_pReal
   forall(n=1:3) Tstar_v(n) = Tstar_v(n) - p_hydro                  ! get deviatoric stress tensor
   
   ! calculate plastic velocity gradient and its tangent according to constitutive law
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
   call constitutive_LpAndItsTangent(Lp_constitutive,dLp_constitutive,Tstar_v,crystallite_Temperature(g,i,e),g,i,e)
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   debug_cumLpCalls = debug_cumLpCalls + 1_pInt
   debug_cumLpTicks  = debug_cumLpTicks + tock-tick
   if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks

   ! update current residuum
   residuum = Lpguess - Lp_constitutive

   ! Check for convergence of loop
   if (.not.(any(residuum/=residuum)) .and. &                       ! exclude any NaN in residuum
       ( maxval(abs(residuum)) < aTol_crystalliteStress .or. &      ! below absolute tolerance .or.
         ( any(abs(crystallite_subdt(g,i,e)*Lpguess) > relevantStrain) .and. & ! worth checking? .and.
             maxval(abs(residuum/Lpguess), abs(crystallite_subdt(g,i,e)*Lpguess) > relevantStrain) < rTol_crystalliteStress & ! below relative tolerance
         ) &
       ) &
      ) &
     exit LpLoop
   
   ! NaN occured at regular speed?
   if (any(residuum/=residuum) .and. leapfrog == 1.0) then
     if (debugger) write(6,*) '::: integrateStress encountered NaN at iteration', NiterationStress
     return

   ! something went wrong at accelerated speed?
   elseif (leapfrog > 1.0_pReal .and. &                             ! at fast pace .and.
            ( sum(residuum*residuum) > sum(residuum_old*residuum_old) .or. &  ! worse residuum .or.
              sum(residuum*residuum_old) < 0.0_pReal .or. &         ! residuum changed sign (overshoot) .or.
              any(residuum/=residuum) &                             ! NaN occured
            ) &
          ) then
     maxleap = 0.5_pReal * leapfrog                                 ! limit next acceleration
     leapfrog = 1.0_pReal                                           ! grinding halt
     
     ! restore old residuum and Lp
     Lpguess = Lpguess_old                                       
     residuum  = residuum_old

   ! residuum got better
   else
     ! calculate Jacobian for correction term 
     dTdLp = 0.0_pReal
     forall (h=1:3,j=1:3,k=1:3,l=1:3,m=1:3,n=1:3) &
       dTdLp(3*(h-1)+j,3*(k-1)+l) = dTdLp(3*(h-1)+j,3*(k-1)+l) + &
                                    C(h,j,l,n)*AB(k,n)+C(h,j,m,l)*BTA(m,k)
     dTdLp = -0.5_pReal*crystallite_subdt(g,i,e)*dTdLp
     dRdLp = math_identity2nd(9) - math_mul99x99(dLp_constitutive,dTdLp)
     invdRdLp = 0.0_pReal
     call math_invert(9,dRdLp,invdRdLp,dummy,error)                 ! invert dR/dLp --> dLp/dR
     if (error) then
       if (debugger) write(6,*) '::: integrateStress failed on dR/dLp inversion at iteration', NiterationStress
       return
     endif

     ! remember current residuum and Lpguess
     residuum_old = residuum
     Lpguess_old = Lpguess 
     
     ! accelerate?
     if (NiterationStress > 1 .and. leapfrog < maxleap) leapfrog = 2.0_pReal * leapfrog
   endif

   ! leapfrog to updated Lp
   forall (k=1:3,l=1:3,m=1:3,n=1:3) & 
     Lpguess(k,l) = Lpguess(k,l) - leapfrog*invdRdLp(3*(k-1)+l,3*(m-1)+n)*residuum(m,n)
 enddo LpLoop

 ! calculate new plastic and elastic deformation gradient
 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new/math_det3x3(invFp_new)**(1.0_pReal/3.0_pReal)  ! regularize by det
 call math_invert3x3(invFp_new,Fp_new,det,error)
 if (error) then
   if (debugger) write(6,*) '::: integrateStress failed on invFp_new inversion at iteration', NiterationStress
   return
 endif
 Fe_new = math_mul33x33(Fg_new,invFp_new)                             ! calc resulting Fe

 ! add volumetric component to 2nd Piola-Kirchhoff stress
 forall (n=1:3) Tstar_v(n) = Tstar_v(n) + p_hydro
 
 ! calculate 1st Piola-Kirchhoff stress
 crystallite_P(:,:,g,i,e) = math_mul33x33(Fe_new,math_mul33x33(math_Mandel6to33(Tstar_v),transpose(invFp_new)))
 
 ! store local values in global variables
 crystallite_Lp(:,:,g,i,e) = Lpguess
 crystallite_Tstar_v(:,g,i,e) = Tstar_v
 crystallite_Fp(:,:,g,i,e) = Fp_new
 crystallite_Fe(:,:,g,i,e) = Fe_new

 ! set return flag to true
 crystallite_integrateStress = .true.
 if (debugger) write(6,*) '::: integrateStress finished at iteration', NiterationStress

!$OMP CRITICAL (distributionStress)
 debug_StressLoopDistribution(NiterationStress) = debug_StressLoopDistribution(NiterationStress) + 1
!$OMP END CRITICAL (distributionStress)

 return

 endfunction
 
 
 
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
 
endfunction


END MODULE
!##############################################################

