
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
 type(p_vec),   dimension(:,:),   allocatable :: homogenization_state0, &        ! pointer array to homogenization state at start of FE increment
                                                 homogenization_subState0, &     ! pointer array to homogenization state at start of homogenization increment
                                                 homogenization_state            ! pointer array to current homogenization state (end of converged time step)
 integer(pInt), dimension(:,:),   allocatable :: homogenization_sizeState, &     ! size of state array per grain
                                                 homogenization_sizePostResults  ! size of postResults array per material point

 real(pReal), dimension(:,:,:,:,:,:), allocatable :: materialpoint_dPdF          ! tangent of first P--K stress at IP
 real(pReal), dimension(:,:,:,:), allocatable :: materialpoint_F0, &             ! def grad of IP at start of FE increment
                                                 materialpoint_F, &              ! def grad of IP to be reached at end of FE increment
                                                 materialpoint_subF0, &          ! def grad of IP at beginning of homogenization increment
                                                 materialpoint_subF, &           ! def grad of IP to be reached at end of homog inc
                                                 materialpoint_P                 ! first P--K stress of IP
 real(pReal), dimension(:,:), allocatable ::     materialpoint_Temperature, &    ! temperature at IP
                                                 materialpoint_subFrac, &
                                                 materialpoint_subStep, &
                                                 materialpoint_subdt

 real(pReal), dimension(:,:,:), allocatable ::   materialpoint_results           ! results array of material point

 logical, dimension(:,:), allocatable ::         materialpoint_requested, &
                                                 materialpoint_converged
 logical, dimension(:,:,:), allocatable ::       materialpoint_doneAndHappy
 integer(pInt) homogenization_maxSizeState,homogenization_maxSizePostResults

CONTAINS

!**************************************
!*      Module initialization         *
!**************************************
subroutine homogenization_init()
 use prec, only: pReal,pInt
 use math, only: math_I3
 use IO, only: IO_error, IO_open_file
 use mesh, only: mesh_maxNips,mesh_NcpElems,mesh_element,FE_Nips
 use material
 use constitutive, only: constitutive_maxSizePostResults
 use crystallite, only: crystallite_Nresults
 use homogenization_isostrain
! use homogenization_RGC

 integer(pInt), parameter :: fileunit = 200
 integer(pInt) e,i,g,myInstance

 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file

 call homogenization_isostrain_init(fileunit)       ! parse all homogenizations of this type

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
 allocate(materialpoint_Temperature(mesh_maxNips,mesh_NcpElems));      materialpoint_Temperature = 0.0_pReal
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
       case (homogenization_isostrain_label)
         if (homogenization_isostrain_sizeState(myInstance) > 0_pInt) then
           allocate(homogenization_state0(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           allocate(homogenization_subState0(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           allocate(homogenization_state(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           homogenization_state0(i,e)%p  = homogenization_isostrain_stateInit(myInstance)
           homogenization_sizeState(i,e) = homogenization_isostrain_sizeState(myInstance)
         endif
         homogenization_sizePostResults(i,e) = homogenization_isostrain_sizePostResults(myInstance)
       case default
         call IO_error(200,ext_msg=homogenization_type(mesh_element(3,e)))      ! unknown type 200 is phase!
     end select
   enddo
 enddo

 homogenization_maxSizeState       = maxval(homogenization_sizeState)
 homogenization_maxSizePostResults = maxval(homogenization_sizePostResults)

 allocate(materialpoint_results(     1+homogenization_maxSizePostResults + &
          homogenization_maxNgrains*(1+crystallite_Nresults+constitutive_maxSizePostResults), mesh_maxNips,mesh_NcpElems))


!    *** Output to MARC output file ***
!$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  homogenization init  -+>>>'
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

end subroutine


!********************************************************************
!*  parallelized calculation of
!*  stress and corresponding tangent
!*  at material points
!********************************************************************
subroutine materialpoint_stressAndItsTangent(&
     updateJaco,&     ! flag to initiate Jacobian updating
     dt &             ! time increment
    )

 use prec,         only: pInt,pReal, subStepMin,nHomog
 use FEsolving,    only: FEsolving_execElem, FEsolving_execIP
 use mesh,         only: mesh_element
 use material,     only: homogenization_Ngrains
 use constitutive, only: constitutive_state0, constitutive_partionedState0, constitutive_state
 use crystallite
 implicit none
 
 real(pReal), intent(in) :: dt
 logical,     intent(in) :: updateJaco
 integer(pInt) homogenization_Niteration
 integer(pInt) g,i,e,myNgrains

! ------ initialize to starting condition ------

 write (6,*)
 write (6,*) 'Material Point start'
 write (6,'(a,/,3(3(f12.7,x)/))') 'F0  of   8 1',materialpoint_F0(1:3,:,8,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'F   of   8 1',materialpoint_F(1:3,:,8,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'Fp0 of 1 8 1',crystallite_Fp0(1:3,:,1,8,1)
 write (6,'(a,/,3(3(f12.7,x)/))') 'Lp0 of 1 8 1',crystallite_Lp0(1:3,:,1,8,1)

!$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
                                                              ! initialize restoration points of grain...
     forall (g = 1:myNgrains) constitutive_partionedState0(g,i,e)%p = constitutive_state0(g,i,e)%p  ! ...microstructures
     crystallite_partionedFp0(:,:,1:myNgrains,i,e)   = crystallite_Fp0(:,:,1:myNgrains,i,e)    ! ...plastic def grads
     crystallite_partionedLp0(:,:,1:myNgrains,i,e)   = crystallite_Lp0(:,:,1:myNgrains,i,e)    ! ...plastic velocity grads
     crystallite_partionedF0(:,:,1:myNgrains,i,e)    = crystallite_F0(:,:,1:myNgrains,i,e)     ! ...def grads
                                                              ! initialize restoration points of ...
     if (homogenization_sizeState(i,e) > 0_pInt) &
       homogenization_subState0(i,e)%p = homogenization_state0(i,e)%p   ! ...internal homogenizaiton state
     materialpoint_subF0(:,:,i,e) = materialpoint_F0(:,:,i,e)           ! ...def grad

     materialpoint_subFrac(i,e) = 0.0_pReal
     materialpoint_subStep(i,e) = 2.0_pReal
     materialpoint_converged(i,e) = .false.                   ! pretend failed step of twice the required size
     materialpoint_requested(i,e) = .true.                    ! everybody requires calculation
   enddo
 enddo
!$OMP END PARALLEL DO


! ------ cutback loop ------

 do while (any(materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMin))  ! cutback loop for material points

!$OMP PARALLEL DO
   do e = FEsolving_execElem(1),FEsolving_execElem(2)           ! iterate over elements to be processed
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)         ! iterate over IPs of this element to be processed
       if (materialpoint_converged(i,e)) then
         materialpoint_subFrac(i,e) = materialpoint_subFrac(i,e) + materialpoint_subStep(i,e)
         materialpoint_subStep(i,e) = min(1.0_pReal-materialpoint_subFrac(i,e), 2.0_pReal * materialpoint_subStep(i,e))
         if (materialpoint_subStep(i,e) > subStepMin) then                       ! still stepping needed
                                                                                 ! wind forward grain starting point of...
           crystallite_partionedF0(:,:,1:myNgrains,i,e)  = crystallite_partionedF(:,:,1:myNgrains,i,e) ! ...def grads
           crystallite_partionedFp0(:,:,1:myNgrains,i,e) = crystallite_Fp(:,:,1:myNgrains,i,e)         ! ...plastic def grads
           crystallite_partionedLp0(:,:,1:myNgrains,i,e) = crystallite_Lp(:,:,1:myNgrains,i,e)         ! ...plastic velocity grads
           forall (g = 1:myNgrains) constitutive_partionedState0(g,i,e)%p = constitutive_state(g,i,e)%p  ! ...microstructures
           if (homogenization_sizeState(i,e) > 0_pInt) &
             homogenization_subState0(i,e)%p = homogenization_state(i,e)%p       ! ...internal state of homog scheme
           materialpoint_subF0(:,:,i,e) = materialpoint_subF(:,:,i,e)            ! ...def grad
         endif
       else
         materialpoint_subStep(i,e) = 0.5_pReal * materialpoint_subStep(i,e)     ! cut step in half and restore...

! ####### why not resetting F0 ?!?!?

         crystallite_Fp(:,:,1:myNgrains,i,e) = crystallite_partionedFp0(:,:,1:myNgrains,i,e)           ! ...plastic def grads
         crystallite_Lp(:,:,1:myNgrains,i,e) = crystallite_partionedLp0(:,:,1:myNgrains,i,e)           ! ...plastic velocity grads
         forall (g = 1:myNgrains) constitutive_state(g,i,e)%p = constitutive_partionedState0(g,i,e)%p  ! ...microstructures
         if (homogenization_sizeState(i,e) > 0_pInt) &
           homogenization_state(i,e)%p = homogenization_subState0(i,e)%p         ! ...internal state of homog scheme
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

! ------ convergence loop material point homogenization ------

   homogenization_Niteration = 0_pInt
   
   do while (any(            materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. materialpoint_doneAndHappy(1,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. homogenization_Niteration < nHomog)           ! convergence loop for materialpoint
     homogenization_Niteration = homogenization_Niteration + 1

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
         endif
       enddo
     enddo
!$OMP END PARALLEL DO

 
     
! --+>> crystallite integration <<+--
!
! based on crystallite_partionedF0,.._partionedF
! incrementing by crystallite_dt

     call crystallite_stressAndItsTangent(updateJaco)                 ! request stress and tangent calculation for constituent grains
     
! --+>> state update <<+--

!$OMP PARALLEL DO
     do e = FEsolving_execElem(1),FEsolving_execElem(2)               ! iterate over elements to be processed
       do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)             ! iterate over IPs of this element to be processed
         if (      materialpoint_requested(i,e) .and. &
             .not. materialpoint_doneAndHappy(1,i,e)) then
           materialpoint_doneAndHappy(:,i,e) = homogenization_updateState(i,e)
           materialpoint_converged(i,e) = all(materialpoint_doneAndHappy(:,i,e))  ! converged if done and happy
         endif
       enddo
     enddo
!$OMP END PARALLEL DO

   enddo                                                           ! homogenization convergence loop  

 enddo                                                             ! cutback loop

 ! check for non-performer: any(.not. converged)
 ! replace with elastic response ?

!$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)                   ! iterate over elements to be processed
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                 ! iterate over IPs of this element to be processed
     call homogenization_averageStressAndItsTangent(i,e)
   enddo
 enddo
!$OMP END PARALLEL DO

 write (6,*) 'Material Point finished'
 write (6,'(a,/,3(3(f12.7,x)/))') 'Lp of 1 8 1',crystallite_Lp(1:3,:,1,8,1)
 
 ! how to deal with stiffness?
 return
 
end subroutine


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
       d = homogenization_sizePostResults(i,e)
       materialpoint_results(c+1,i,e) = d; c = c+1_pInt         ! tell size of homogenization results
       if (d > 0_pInt) then                                        ! any homogenization results to mention?
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

 end subroutine
 
 
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

 implicit none
 
 integer(pInt), intent(in) :: ip,el
 
 select case(homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label)
     call homogenization_isostrain_partitionDeformation(crystallite_partionedF(:,:,:,ip,el), &
                                                        crystallite_partionedF0(:,:,:,ip,el),&
                                                        materialpoint_subF(:,:,ip,el),&
                                                        homogenization_state(ip,el),ip,el)
 end select

end subroutine


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
 use crystallite, only: crystallite_P,crystallite_dPdF

 use homogenization_isostrain
 implicit none
 
 integer(pInt), intent(in) :: ip,el
 logical, dimension(2) :: homogenization_updateState
 
 select case(homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label)
     homogenization_updateState = &
     homogenization_isostrain_updateState(homogenization_state(ip,el), &
                                          crystallite_P(:,:,:,ip,el),crystallite_dPdF(:,:,:,:,:,ip,el),ip,el)
 end select

 return
 
end function


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

 use homogenization_isostrain
 implicit none
 
 integer(pInt), intent(in) :: ip,el
 
 select case(homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label)
     call homogenization_isostrain_averageStressAndItsTangent(materialpoint_P(:,:,ip,el), materialpoint_dPdF(:,:,:,:,ip,el),&
                                                              crystallite_P(:,:,:,ip,el),crystallite_dPdF(:,:,:,:,:,ip,el),ip,el)
 end select

 return
 
end subroutine


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
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ip,el
 real(pReal), dimension(homogenization_sizePostResults(ip,el)) :: homogenization_postResults

 homogenization_postResults = 0.0_pReal
 select case (homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label)
     homogenization_postResults = homogenization_isostrain_postResults(homogenization_state(ip,el),ip,el)

 end select

 return

end function

END MODULE