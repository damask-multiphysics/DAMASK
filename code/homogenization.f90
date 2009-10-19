!* $Id$
!***************************************
!*      Module: HOMOGENIZATION         *
!***************************************
!* contains:                           *
!* - _init                             *
!* - materialpoint_stressAndItsTangent *
!* - _partitionDeformation             *
!* - _updateState                      *
!* - _averageStressAndItsTangent       *
!* - _postResults                      *
!***************************************

MODULE homogenization

!*** Include other modules ***
 use prec, only: pInt,pReal,p_vec
 implicit none

! ****************************************************************
! *** General variables for the homogenization at a            ***
! *** material point                                           ***
! ****************************************************************
 type(p_vec), dimension(:,:), allocatable ::         homogenization_state0, &        ! pointer array to homogenization state at start of FE increment
                                                     homogenization_subState0, &     ! pointer array to homogenization state at start of homogenization increment
                                                     homogenization_state            ! pointer array to current homogenization state (end of converged time step)
 integer(pInt), dimension(:,:), allocatable ::       homogenization_sizeState, &     ! size of state array per grain
                                                     homogenization_sizePostResults  ! size of postResults array per material point

 real(pReal), dimension(:,:,:,:,:,:), allocatable :: materialpoint_dPdF              ! tangent of first P--K stress at IP
 real(pReal), dimension(:,:,:,:), allocatable ::     materialpoint_F0, &             ! def grad of IP at start of FE increment
                                                     materialpoint_F, &              ! def grad of IP to be reached at end of FE increment
                                                     materialpoint_subF0, &          ! def grad of IP at beginning of homogenization increment
                                                     materialpoint_subF, &           ! def grad of IP to be reached at end of homog inc
                                                     materialpoint_P                 ! first P--K stress of IP
 real(pReal), dimension(:,:), allocatable ::         materialpoint_Temperature, &    ! temperature at IP
                                                     materialpoint_subFrac, &
                                                     materialpoint_subStep, &
                                                     materialpoint_subdt

 real(pReal), dimension(:,:,:), allocatable ::       materialpoint_results           ! results array of material point

 logical, dimension(:,:), allocatable ::             materialpoint_requested, &
                                                     materialpoint_converged
 logical, dimension(:,:,:), allocatable ::           materialpoint_doneAndHappy
 integer(pInt)                                       homogenization_maxSizeState, &
                                                     homogenization_maxSizePostResults, &
                                                     materialpoint_sizeResults

CONTAINS

!**************************************
!*      Module initialization         *
!**************************************
subroutine homogenization_init(Temperature)
 use prec, only: pReal,pInt
 use math, only: math_I3
 use IO, only: IO_error, IO_open_file
 use mesh, only: mesh_maxNips,mesh_NcpElems,mesh_element,FE_Nips
 use material
 use constitutive, only: constitutive_maxSizePostResults
 use crystallite, only: crystallite_Nresults
 use homogenization_isostrain
 use homogenization_RGC                             ! RGC homogenization added <<<updated 31.07.2009>>>

 real(pReal) Temperature
 integer(pInt), parameter :: fileunit = 200
 integer(pInt) e,i,g,myInstance,j

 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file

 call homogenization_isostrain_init(fileunit)       ! parse all homogenizations of this type
 call homogenization_RGC_init(fileunit)             ! RGC homogenization added <<<updated 31.07.2009>>>

 close(fileunit)

 allocate(homogenization_state0(mesh_maxNips,mesh_NcpElems))
 allocate(homogenization_subState0(mesh_maxNips,mesh_NcpElems))
 allocate(homogenization_state(mesh_maxNips,mesh_NcpElems))
 allocate(homogenization_sizeState(mesh_maxNips,mesh_NcpElems));       homogenization_sizeState = 0_pInt
 allocate(homogenization_sizePostResults(mesh_maxNips,mesh_NcpElems)); homogenization_sizePostResults = 0_pInt

 allocate(materialpoint_dPdF(3,3,3,3,mesh_maxNips,mesh_NcpElems));     materialpoint_dPdF    = 0.0_pReal
 allocate(materialpoint_F0(3,3,mesh_maxNips,mesh_NcpElems));
 allocate(materialpoint_F(3,3,mesh_maxNips,mesh_NcpElems));            materialpoint_F       = 0.0_pReal
 allocate(materialpoint_subF0(3,3,mesh_maxNips,mesh_NcpElems));        materialpoint_subF0   = 0.0_pReal
 allocate(materialpoint_subF(3,3,mesh_maxNips,mesh_NcpElems));         materialpoint_subF    = 0.0_pReal
 allocate(materialpoint_P(3,3,mesh_maxNips,mesh_NcpElems));            materialpoint_P       = 0.0_pReal
 allocate(materialpoint_Temperature(mesh_maxNips,mesh_NcpElems));      materialpoint_Temperature = Temperature
 allocate(materialpoint_subFrac(mesh_maxNips,mesh_NcpElems));          materialpoint_subFrac = 0.0_pReal
 allocate(materialpoint_subStep(mesh_maxNips,mesh_NcpElems));          materialpoint_subStep = 0.0_pReal
 allocate(materialpoint_subdt(mesh_maxNips,mesh_NcpElems));            materialpoint_subdt   = 0.0_pReal
 allocate(materialpoint_requested(mesh_maxNips,mesh_NcpElems));        materialpoint_requested = .false.
 allocate(materialpoint_converged(mesh_maxNips,mesh_NcpElems));        materialpoint_converged = .true.
 allocate(materialpoint_doneAndHappy(2,mesh_maxNips,mesh_NcpElems));   materialpoint_doneAndHappy = .true.

 forall (i = 1:mesh_maxNips,e = 1:mesh_NcpElems)
   materialpoint_F0(:,:,i,e) = math_I3
   materialpoint_F(:,:,i,e)  = math_I3
 end forall

 do e = 1,mesh_NcpElems                                  ! loop over elements
   myInstance = homogenization_typeInstance(mesh_element(3,e))
   do i = 1,FE_Nips(mesh_element(2,e))                   ! loop over IPs
     select case(homogenization_type(mesh_element(3,e)))
!* isostrain    
       case (homogenization_isostrain_label)
         if (homogenization_isostrain_sizeState(myInstance) > 0_pInt) then
           allocate(homogenization_state0(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           allocate(homogenization_subState0(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           allocate(homogenization_state(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           homogenization_state0(i,e)%p  = homogenization_isostrain_stateInit(myInstance)
           homogenization_sizeState(i,e) = homogenization_isostrain_sizeState(myInstance)
         endif
         homogenization_sizePostResults(i,e) = homogenization_isostrain_sizePostResults(myInstance)
!* RGC homogenization: added <<<updated 31.07.2009>>>
       case (homogenization_RGC_label)
         if (homogenization_RGC_sizeState(myInstance) > 0_pInt) then
           allocate(homogenization_state0(i,e)%p(homogenization_RGC_sizeState(myInstance)))
           allocate(homogenization_subState0(i,e)%p(homogenization_RGC_sizeState(myInstance)))
           allocate(homogenization_state(i,e)%p(homogenization_RGC_sizeState(myInstance)))
           homogenization_state0(i,e)%p  = homogenization_RGC_stateInit(myInstance)
           homogenization_sizeState(i,e) = homogenization_RGC_sizeState(myInstance)
         endif
         homogenization_sizePostResults(i,e) = homogenization_RGC_sizePostResults(myInstance)
       case default
         call IO_error(201,ext_msg=homogenization_type(mesh_element(3,e)))      ! unknown type 201 is homogenization!
     end select
   enddo
 enddo

 homogenization_maxSizeState       = maxval(homogenization_sizeState)
 homogenization_maxSizePostResults = maxval(homogenization_sizePostResults)

 materialpoint_sizeResults = 1+ 1+homogenization_maxSizePostResults + &    ! grain count, homogSize, homogResult
          homogenization_maxNgrains*(1+crystallite_Nresults+constitutive_maxSizePostResults)
 allocate(materialpoint_results(  materialpoint_sizeResults, mesh_maxNips,mesh_NcpElems))


!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  homogenization init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'homogenization_state0:          ', shape(homogenization_state0)
 write(6,'(a32,x,7(i5,x))') 'homogenization_subState0:       ', shape(homogenization_subState0)
 write(6,'(a32,x,7(i5,x))') 'homogenization_state:           ', shape(homogenization_state)
 write(6,'(a32,x,7(i5,x))') 'homogenization_sizeState:       ', shape(homogenization_sizeState)
 write(6,'(a32,x,7(i5,x))') 'homogenization_sizePostResults: ', shape(homogenization_sizePostResults)
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_dPdF:             ', shape(materialpoint_dPdF)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_F0:               ', shape(materialpoint_F0)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_F:                ', shape(materialpoint_F)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_subF0:            ', shape(materialpoint_subF0)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_subF:             ', shape(materialpoint_subF)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_P:                ', shape(materialpoint_P)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_Temperature:      ', shape(materialpoint_Temperature)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_subFrac:          ', shape(materialpoint_subFrac)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_subStep:          ', shape(materialpoint_subStep)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_subdt:            ', shape(materialpoint_subdt)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_requested:        ', shape(materialpoint_requested)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_converged:        ', shape(materialpoint_converged)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_doneAndHappy:     ', shape(materialpoint_doneAndHappy)
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'materialpoint_results:          ', shape(materialpoint_results)
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'maxSizeState:       ', homogenization_maxSizeState
 write(6,'(a32,x,7(i5,x))') 'maxSizePostResults: ', homogenization_maxSizePostResults
 call flush(6)
!$OMP END CRITICAL (write2out)

 return

endsubroutine


!********************************************************************
!*  parallelized calculation of
!*  stress and corresponding tangent
!*  at material points
!********************************************************************
subroutine materialpoint_stressAndItsTangent(&
     updateJaco,&     ! flag to initiate Jacobian updating
     dt &             ! time increment
    )

 use prec, only:          pInt, &
                          pReal
 use numerics, only:      subStepMin, &
                          nHomog, &
                          nMPstate
 use FEsolving, only:     FEsolving_execElem, &
                          FEsolving_execIP, &
                          terminallyIll
 use mesh, only:          mesh_element
 use material, only:      homogenization_Ngrains
 use constitutive, only:  constitutive_state0, &
                          constitutive_partionedState0, &
                          constitutive_state
 use crystallite, only:   crystallite_Temperature, &
                          crystallite_F0, &
                          crystallite_Fp0, &
                          crystallite_Fp, &
                          crystallite_Lp0, &
                          crystallite_Lp, &
                          crystallite_Tstar0_v, &
                          crystallite_Tstar_v, &
                          crystallite_partionedTemperature0, &
                          crystallite_partionedF0, &
                          crystallite_partionedF, &
                          crystallite_partionedFp0, &
                          crystallite_partionedLp0, &
                          crystallite_partionedTstar0_v, &
                          crystallite_dt, &
                          crystallite_requested, &
                          crystallite_converged, &
                          crystallite_stressAndItsTangent
 use debug, only:         debugger, &
                          debug_MaterialpointLoopDistribution, &
                          debug_MaterialpointStateLoopDistribution
                          
 implicit none
 
 real(pReal), intent(in) :: dt
 logical,     intent(in) :: updateJaco
 integer(pInt) NiterationHomog,NiterationMPstate
 integer(pInt) g,i,e,myNgrains

! ------ initialize to starting condition ------

 write (6,*)
 write (6,*) 'Material Point start'
 write (6,'(a,/,(f12.7,x))')      'Temp0  of   1 1'  ,materialpoint_Temperature(1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'F0     of   1 1',materialpoint_F0(1:3,:,1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'F      of   1 1',materialpoint_F(1:3,:,1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'Fp0    of 1 1 1',crystallite_Fp0(1:3,:,1,1,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'Lp0    of 1 1 1',crystallite_Lp0(1:3,:,1,1,1)

!$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                 ! iterate over elements to be processed
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                               ! iterate over IPs of this element to be processed
     
     ! initialize restoration points of grain...
     forall (g = 1:myNgrains) constitutive_partionedState0(g,i,e)%p = constitutive_state0(g,i,e)%p   ! ...microstructures
     crystallite_partionedTemperature0(1:myNgrains,i,e) = materialpoint_Temperature(i,e)             ! ...temperatures
     crystallite_partionedFp0(:,:,1:myNgrains,i,e)   = crystallite_Fp0(:,:,1:myNgrains,i,e)          ! ...plastic def grads
     crystallite_partionedLp0(:,:,1:myNgrains,i,e)   = crystallite_Lp0(:,:,1:myNgrains,i,e)          ! ...plastic velocity grads
     crystallite_partionedF0(:,:,1:myNgrains,i,e)    = crystallite_F0(:,:,1:myNgrains,i,e)           ! ...def grads
     crystallite_partionedTstar0_v(:,1:myNgrains,i,e)= crystallite_Tstar0_v(:,1:myNgrains,i,e)       ! ...2nd PK stress
     
     ! initialize restoration points of ...
     if (homogenization_sizeState(i,e) > 0_pInt) &
       homogenization_subState0(i,e)%p = homogenization_state0(i,e)%p                               ! ...internal homogenization state
     materialpoint_subF0(:,:,i,e) = materialpoint_F0(:,:,i,e)                                       ! ...def grad

     materialpoint_subFrac(i,e) = 0.0_pReal
     materialpoint_subStep(i,e) = 8.0_pReal
     materialpoint_converged(i,e) = .false.                                                         ! pretend failed step of twice the required size
     materialpoint_requested(i,e) = .true.                                                          ! everybody requires calculation
   enddo
 enddo
!$OMP END PARALLEL DO

 NiterationHomog = 0_pInt
 
! ------ cutback loop ------

 do while (any(materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMin))  ! cutback loop for material points

!    write(6,'(a,/,125(8(f8.5,x),/))') 'mp_subSteps',materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2))
!$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)                                               ! iterate over elements to be processed
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                             ! iterate over IPs of this element to be processed
       
       debugger = (e == 1 .and. i == 1)
       
       ! if our materialpoint converged then we are either finished or have to wind forward
       if (materialpoint_converged(i,e)) then
         if (debugger) then
           !$OMP CRITICAL (write2out)
             write(6,'(a21,f10.8,a34,f10.8,a37,/)') 'winding forward from ', &
               materialpoint_subFrac(i,e), ' to current materialpoint_subFrac ', &
               materialpoint_subFrac(i,e)+materialpoint_subStep(i,e),' in materialpoint_stressAndItsTangent'
           !$OMPEND CRITICAL (write2out)
         endif
         
         ! calculate new subStep and new subFrac
         materialpoint_subFrac(i,e) = materialpoint_subFrac(i,e) + materialpoint_subStep(i,e)
         materialpoint_subStep(i,e) = min(1.0_pReal-materialpoint_subFrac(i,e), 1.0_pReal * materialpoint_subStep(i,e))   ! keep cut back time step (no acceleration)
                  
         ! still stepping needed
         if (materialpoint_subStep(i,e) > subStepMin) then
         
           ! wind forward grain starting point of...
           crystallite_partionedTemperature0(1:myNgrains,i,e) = crystallite_Temperature(1:myNgrains,i,e)  ! ...temperatures
           crystallite_partionedF0(:,:,1:myNgrains,i,e)     = crystallite_partionedF(:,:,1:myNgrains,i,e) ! ...def grads
           crystallite_partionedFp0(:,:,1:myNgrains,i,e)    = crystallite_Fp(:,:,1:myNgrains,i,e)         ! ...plastic def grads
           crystallite_partionedLp0(:,:,1:myNgrains,i,e)    = crystallite_Lp(:,:,1:myNgrains,i,e)         ! ...plastic velocity grads
           crystallite_partionedTstar0_v(:,1:myNgrains,i,e) = crystallite_Tstar_v(:,1:myNgrains,i,e)      ! ...2nd PK stress
           forall (g = 1:myNgrains) constitutive_partionedState0(g,i,e)%p = constitutive_state(g,i,e)%p   ! ...microstructures
           if (homogenization_sizeState(i,e) > 0_pInt) &
             homogenization_subState0(i,e)%p = homogenization_state(i,e)%p                                ! ...internal state of homog scheme
           materialpoint_subF0(:,:,i,e) = materialpoint_subF(:,:,i,e)                                     ! ...def grad
         elseif (materialpoint_requested(i,e)) then                                                       ! this materialpoint just converged    ! already at final time (??)
          !$OMP CRITICAL (distributionHomog)
            debug_MaterialpointLoopDistribution(min(nHomog+1,NiterationHomog)) = &
              debug_MaterialpointLoopDistribution(min(nHomog+1,NiterationHomog)) + 1
          !$OMPEND CRITICAL (distributionHomog)
         endif
       
       ! materialpoint didn't converge, so we need a cutback here
       else
       
         materialpoint_subStep(i,e) = 0.125_pReal * materialpoint_subStep(i,e)                            ! crystallite had severe trouble, so do a significant cutback
         
         if (debugger) then
           !$OMP CRITICAL (write2out)
             write(6,'(a82,f10.8,/)') 'cutback step in materialpoint_stressAndItsTangent with new materialpoint_subStep: ',&
                                       materialpoint_subStep(i,e)
           !$OMPEND CRITICAL (write2out)
         endif

         ! restore...
         crystallite_Temperature(1:myNgrains,i,e) = crystallite_partionedTemperature0(1:myNgrains,i,e)    ! ...temperatures
         crystallite_Fp(:,:,1:myNgrains,i,e)    = crystallite_partionedFp0(:,:,1:myNgrains,i,e)           ! ...plastic def grads
         crystallite_Lp(:,:,1:myNgrains,i,e)    = crystallite_partionedLp0(:,:,1:myNgrains,i,e)           ! ...plastic velocity grads
         crystallite_Tstar_v(:,1:myNgrains,i,e) = crystallite_partionedTstar0_v(:,1:myNgrains,i,e)        ! ...2nd PK stress
         forall (g = 1:myNgrains) constitutive_state(g,i,e)%p = constitutive_partionedState0(g,i,e)%p     ! ...microstructures
         if (homogenization_sizeState(i,e) > 0_pInt) &
           homogenization_state(i,e)%p = homogenization_subState0(i,e)%p                                  ! ...internal state of homog scheme
       
       endif
     
       materialpoint_requested(i,e) = materialpoint_subStep(i,e) > subStepMin
       if (materialpoint_requested(i,e)) then
         materialpoint_subF(:,:,i,e) = materialpoint_subF0(:,:,i,e) + &
                                       materialpoint_subStep(i,e) * (materialpoint_F(:,:,i,e) - materialpoint_F0(:,:,i,e))
         materialpoint_subdt(i,e)    = materialpoint_subStep(i,e) * dt
         materialpoint_doneAndHappy(:,i,e) = (/.false.,.true./)
       endif
     enddo
   enddo
!$OMP END PARALLEL DO

!* Checks for cutback/substepping loops: added <<<updated 31.07.2009>>>
 ! write (6,'(a,/,8(L,x))') 'MP exceeds substep min',materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMin
 ! write (6,'(a,/,8(L,x))') 'MP requested',materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2))
 ! write (6,'(a,/,8(f6.4,x))') 'MP subFrac',materialpoint_subFrac(:,FEsolving_execELem(1):FEsolving_execElem(2))
 ! write (6,'(a,/,8(f6.4,x))') 'MP subStep',materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2))

! ------ convergence loop material point homogenization ------

   NiterationMPstate = 0_pInt
   
   do while (any(            materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. materialpoint_doneAndHappy(1,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. NiterationMPstate < nMPstate)           ! convergence loop for materialpoint
     NiterationMPstate = NiterationMPstate + 1

!      write(6,'(a,/,125(8(l,x),/))') 'material point request and not done', &
!                                     materialpoint_requested .and. .not. materialpoint_doneAndHappy(1,:,:)

! --+>> deformation partitioning <<+--
!
! based on materialpoint_subF0,.._subF,
!          crystallite_partionedF0,
!          homogenization_state
! results in crystallite_partionedF

!$OMP PARALLEL DO
     do e = FEsolving_execElem(1),FEsolving_execElem(2)               ! iterate over elements to be processed
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)             ! iterate over IPs of this element to be processed
         if (      materialpoint_requested(i,e) .and. &               ! process requested but...
             .not. materialpoint_doneAndHappy(1,i,e)) then            ! ...not yet done material points
           call homogenization_partitionDeformation(i,e)              ! partition deformation onto constituents
           crystallite_dt(1:myNgrains,i,e) = materialpoint_subdt(i,e) ! propagate materialpoint dt to grains
           crystallite_requested(1:myNgrains,i,e) = .true.            ! request calculation for constituents
         else
           crystallite_requested(1:myNgrains,i,e) = .false.           ! calculation for constituents not required anymore
         endif
       enddo
     enddo
!$OMP END PARALLEL DO
!      write(6,'(a,/,125(8(8(l,x),2x),/))') 'crystallite request with updated partitioning', crystallite_requested
 
     
! --+>> crystallite integration <<+--
!
! based on crystallite_partionedF0,.._partionedF
! incrementing by crystallite_dt
     call crystallite_stressAndItsTangent(updateJaco)                 ! request stress and tangent calculation for constituent grains

!      write(6,'(a,/,125(8(8(l,x),2x),/))') 'crystallite converged', crystallite_converged
     
! --+>> state update <<+--

!$OMP PARALLEL DO
     do e = FEsolving_execElem(1),FEsolving_execElem(2)               ! iterate over elements to be processed
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)             ! iterate over IPs of this element to be processed
         if (      materialpoint_requested(i,e) .and. &
             .not. materialpoint_doneAndHappy(1,i,e)) then
           if (.not. all(crystallite_converged(:,i,e))) then
             materialpoint_doneAndHappy(:,i,e) = (/.true.,.false./)
           else
             materialpoint_doneAndHappy(:,i,e) = homogenization_updateState(i,e)
           endif
           materialpoint_converged(i,e) = all(materialpoint_doneAndHappy(:,i,e))  ! converged if done and happy
           if (materialpoint_converged(i,e)) &                                    ! added <<<updated 31.07.2009>>>
             debug_MaterialpointStateLoopdistribution(NiterationMPstate) = &
             debug_MaterialpointStateLoopdistribution(NiterationMPstate) + 1
         endif
       enddo
     enddo
!$OMP END PARALLEL DO
!      write(6,'(a,/,125(8(l,x),/))') 'material point done', materialpoint_doneAndHappy(1,:,:)
!      write(6,'(a,/,125(8(l,x),/))') 'material point converged', materialpoint_converged

   enddo                                                           ! homogenization convergence loop  

   NiterationHomog = NiterationHomog +1_pInt

 enddo                                                             ! cutback loop

 ! check for non-performer: any(.not. converged)
 ! replace everybody with odd response ?

!$OMP PARALLEL DO
elementLoop: do e = FEsolving_execElem(1),FEsolving_execElem(2)       ! iterate over elements to be processed
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                 ! iterate over IPs of this element to be processed
     if (materialpoint_converged(i,e)) then
       call homogenization_averageStressAndItsTangent(i,e)
       call homogenization_averageTemperature(i,e)	 
     else
       terminallyIll = .true.
       exit elementLoop
     endif
   enddo
 enddo elementLoop
!$OMP END PARALLEL DO

 write (6,*)
 write (6,*) 'Material Point end'
 write (6,*)
 return
 
endsubroutine


!********************************************************************
!*  parallelized calculation of
!*  result array at material points
!********************************************************************
subroutine materialpoint_postResults(dt)

 use FEsolving,    only: FEsolving_execElem, FEsolving_execIP
 use mesh,         only: mesh_element
 use material,     only: homogenization_Ngrains
 use constitutive, only: constitutive_sizePostResults, constitutive_postResults
 use crystallite
 implicit none

 real(pReal), intent(in) :: dt
 integer(pInt) g,i,e,c,d,myNgrains

!$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
       c = 0_pInt
       materialpoint_results(c+1,i,e) = myNgrains; c = c+1_pInt ! tell number of grains at materialpoint
       d = homogenization_sizePostResults(i,e)
       materialpoint_results(c+1,i,e) = d; c = c+1_pInt         ! tell size of homogenization results
       if (d > 0_pInt) then                                     ! any homogenization results to mention?
         materialpoint_results(c+1:c+d,i,e) = &                 ! tell homogenization results
         homogenization_postResults(i,e);  c = c+d
       endif
       do g = 1,myNgrains                                       ! 
         d = crystallite_Nresults+constitutive_sizePostResults(g,i,e)
         materialpoint_results(c+1,i,e) = d; c = c+1_pInt       ! tell size of crystallite results
         materialpoint_results(c+1:c+d,i,e) = &                 ! tell crystallite results
           crystallite_postResults(crystallite_Tstar_v(:,g,i,e),crystallite_Temperature(g,i,e),dt,g,i,e); c = c+d
       enddo
     enddo
   enddo
!$OMP END PARALLEL DO

 endsubroutine
 
 
!********************************************************************
! partition material point def grad onto constituents
!********************************************************************
subroutine homogenization_partitionDeformation(&
   ip, &            ! integration point
   el  &            ! element
  )

 use prec,        only: pReal,pInt
 use mesh,        only: mesh_element
 use material,    only: homogenization_type, homogenization_maxNgrains
 use crystallite, only: crystallite_partionedF0,crystallite_partionedF
 use homogenization_isostrain
 use homogenization_RGC                          ! RGC homogenization added <<<updated 31.07.2009>>>

 implicit none
 
 integer(pInt), intent(in) :: ip,el
 
 select case(homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label)
!* isostrain
     call homogenization_isostrain_partitionDeformation(crystallite_partionedF(:,:,:,ip,el), &
                                                        crystallite_partionedF0(:,:,:,ip,el),&
                                                        materialpoint_subF(:,:,ip,el),&
                                                        homogenization_state(ip,el), &
                                                        ip, &
                                                        el)
!* RGC homogenization added <<<updated 31.07.2009>>>
   case (homogenization_RGC_label)
     call homogenization_RGC_partitionDeformation(crystallite_partionedF(:,:,:,ip,el), &
                                                  crystallite_partionedF0(:,:,:,ip,el),&
                                                  materialpoint_subF(:,:,ip,el),&
                                                  homogenization_state(ip,el), &
                                                  ip, &
                                                  el)
 end select

endsubroutine


!********************************************************************
! update the internal state of the homogenization scheme
! and tell whether "done" and "happy" with result
!********************************************************************
function homogenization_updateState(&
   ip, &            ! integration point
   el  &            ! element
  )
 use prec,        only: pReal,pInt
 use mesh,        only: mesh_element
 use material,    only: homogenization_type, homogenization_maxNgrains
 use crystallite, only: crystallite_P,crystallite_dPdF,crystallite_partionedF,crystallite_partionedF0  ! modified <<<updated 31.07.2009>>>

 use homogenization_isostrain
 use homogenization_RGC              ! RGC homogenization added <<<updated 31.07.2009>>>
 implicit none
 
 integer(pInt), intent(in) :: ip,el
 logical, dimension(2) :: homogenization_updateState
 
 select case(homogenization_type(mesh_element(3,el)))
!* isostrain
   case (homogenization_isostrain_label)
     homogenization_updateState = homogenization_isostrain_updateState( homogenization_state(ip,el), &
                                                                        crystallite_P(:,:,:,ip,el), &
                                                                        crystallite_dPdF(:,:,:,:,:,ip,el), &
                                                                        ip, &
                                                                        el)
!* RGC homogenization added <<<updated 31.07.2009>>>
   case (homogenization_RGC_label)
     homogenization_updateState = homogenization_RGC_updateState( homogenization_state(ip,el), &
                                                                  crystallite_P(:,:,:,ip,el), &
                                                                  crystallite_partionedF(:,:,:,ip,el), &
                                                                  crystallite_partionedF0(:,:,:,ip,el),&
                                                                  materialpoint_subF(:,:,ip,el),&
                                                                  crystallite_dPdF(:,:,:,:,:,ip,el), &
                                                                  ip, &
                                                                  el)
 end select

 return
 
endfunction


!********************************************************************
! derive average stress and stiffness from constituent quantities
!********************************************************************
subroutine homogenization_averageStressAndItsTangent(&
   ip, &            ! integration point
   el  &            ! element
  )
 use prec,        only: pReal,pInt
 use mesh,        only: mesh_element
 use material,    only: homogenization_type, homogenization_maxNgrains
 use crystallite, only: crystallite_P,crystallite_dPdF

 use homogenization_RGC            ! RGC homogenization added <<<updated 31.07.2009>>>
 use homogenization_isostrain
 implicit none
 
 integer(pInt), intent(in) :: ip,el
 
 select case(homogenization_type(mesh_element(3,el)))
!* isostrain
   case (homogenization_isostrain_label)
     call homogenization_isostrain_averageStressAndItsTangent( materialpoint_P(:,:,ip,el), &
                                                               materialpoint_dPdF(:,:,:,:,ip,el),&
                                                               crystallite_P(:,:,:,ip,el), &
                                                               crystallite_dPdF(:,:,:,:,:,ip,el), &
                                                               ip, &
                                                               el)
!* RGC homogenization added <<<updated 31.07.2009>>>
   case (homogenization_RGC_label)
     call homogenization_RGC_averageStressAndItsTangent( materialpoint_P(:,:,ip,el), &
                                                         materialpoint_dPdF(:,:,:,:,ip,el),&
                                                         crystallite_P(:,:,:,ip,el), &
                                                         crystallite_dPdF(:,:,:,:,:,ip,el), &
                                                         ip, &
                                                         el)
 end select

 return
 
endsubroutine


!********************************************************************
! derive average stress and stiffness from constituent quantities
!********************************************************************
subroutine homogenization_averageTemperature(&
   ip, &            ! integration point
   el  &            ! element
  )
 use prec,        only: pReal,pInt
 use mesh,        only: mesh_element
 use material,    only: homogenization_type, homogenization_maxNgrains
 use crystallite, only: crystallite_Temperature

 use homogenization_isostrain
 use homogenization_RGC               ! RGC homogenization added <<<updated 31.07.2009>>>
 implicit none
 
 integer(pInt), intent(in) :: ip,el
 
 select case(homogenization_type(mesh_element(3,el)))
!* isostrain
   case (homogenization_isostrain_label)
     materialpoint_Temperature(ip,el) =  homogenization_isostrain_averageTemperature(crystallite_Temperature(:,ip,el), ip, el)
!* RGC homogenization added <<<updated 31.07.2009>>>
   case (homogenization_RGC_label)
     materialpoint_Temperature(ip,el) =  homogenization_RGC_averageTemperature(crystallite_Temperature(:,ip,el), ip, el)
 end select

 return
 
endsubroutine


!********************************************************************
! return array of homogenization results for post file inclusion
! call only, if homogenization_sizePostResults(ip,el) > 0 !!
!********************************************************************
function homogenization_postResults(&
   ip, &            ! integration point
   el  &            ! element
  )
 use prec,     only: pReal,pInt
 use mesh,     only: mesh_element
 use material, only: homogenization_type
 use homogenization_isostrain
 use homogenization_RGC             ! RGC homogenization added <<<updated 31.07.2009>>>
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ip,el
 real(pReal), dimension(homogenization_sizePostResults(ip,el)) :: homogenization_postResults

 homogenization_postResults = 0.0_pReal
 select case (homogenization_type(mesh_element(3,el)))
!* isostrain
   case (homogenization_isostrain_label)
     homogenization_postResults = homogenization_isostrain_postResults(homogenization_state(ip,el),ip,el)
!* RGC homogenization added <<<updated 31.07.2009>>>
   case (homogenization_RGC_label)
     homogenization_postResults = homogenization_RGC_postResults(homogenization_state(ip,el),ip,el)
 end select

 return

endfunction

END MODULE
