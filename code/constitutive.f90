!* $Id$
!************************************
!*      Module: CONSTITUTIVE        *
!************************************
!* contains:                        *
!* - constitutive equations         *
!* - parameters definition          *
!************************************

MODULE constitutive

!*** Include other modules ***
 use prec
 implicit none

 type(p_vec),   dimension(:,:,:), allocatable :: constitutive_state0, &          ! pointer array to microstructure at start of FE inc
                                                 constitutive_partionedState0, & ! pointer array to microstructure at start of homogenization inc
                                                 constitutive_subState0, &       ! pointer array to microstructure at start of crystallite inc
                                                 constitutive_state, &           ! pointer array to current microstructure (end of converged time step)
                                                 constitutive_dotState, &        ! pointer array to evolution of current microstructure
                                                 constitutive_previousDotState, &! pointer array to previous evolution of current microstructure
                                                 constitutive_previousDotState2, &! pointer array to previous evolution of current microstructure
                                                 constitutive_relevantState      ! relevant state values
 integer(pInt), dimension(:,:,:), allocatable :: constitutive_sizeDotState, &    ! size of dotState array
                                                 constitutive_sizeState, &       ! size of state array per grain
                                                 constitutive_sizePostResults    ! size of postResults array per grain
 integer(pInt)                                   constitutive_maxSizeDotState, &
                                                 constitutive_maxSizeState, &
                                                 constitutive_maxSizePostResults

CONTAINS
!****************************************
!* - constitutive_init
!* - constitutive_homogenizedC
!* - constitutive_averageBurgers
!* - constitutive_microstructure
!* - constitutive_LpAndItsTangent
!* - constitutive_collectDotState
!* - constitutive_collectDotTemperature
!* - constitutive_postResults
!****************************************


subroutine constitutive_init()
!**************************************
!*      Module initialization         *
!**************************************
 use prec, only: pReal,pInt
 use debug, only: debugger, selectiveDebugger, debug_e, debug_i, debug_g
 use IO, only: IO_error, IO_open_file, IO_open_jobFile
 use mesh, only: mesh_maxNips,mesh_NcpElems,mesh_element,FE_Nips
 use material
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal

 integer(pInt), parameter :: fileunit = 200
 integer(pInt) e,i,g,p,myInstance,myNgrains
 integer(pInt), dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 logical knownConstitution
 
 if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file

 call constitutive_j2_init(fileunit)                     ! parse all phases of this constitution
 call constitutive_phenopowerlaw_init(fileunit)
 call constitutive_dislotwin_init(fileunit)
 call constitutive_nonlocal_init(fileunit)


 close(fileunit)

! write description file for constitutive phase output

 if(.not. IO_open_jobFile(fileunit,'outputConstitutive')) call IO_error (50) ! problems in writing file
 
 do p = 1,material_Nphase
   i = phase_constitutionInstance(p)                     ! which instance of a constitution is present phase
   knownConstitution = .true.                            ! assume valid
   select case(phase_constitution(p))                    ! split per constitiution
     case (constitutive_j2_label)
       thisOutput => constitutive_j2_output
       thisSize   => constitutive_j2_sizePostResult
     case (constitutive_phenopowerlaw_label)
       thisOutput => constitutive_phenopowerlaw_output
       thisSize   => constitutive_phenopowerlaw_sizePostResult
     case (constitutive_dislotwin_label)
       thisOutput => constitutive_dislotwin_output
       thisSize   => constitutive_dislotwin_sizePostResult
     case (constitutive_nonlocal_label)
       thisOutput => constitutive_nonlocal_output
       thisSize   => constitutive_nonlocal_sizePostResult
     case default
       knownConstitution = .false.
   end select   

   write(fileunit,*)
   write(fileunit,'(a)') '['//trim(phase_name(p))//']'
   write(fileunit,*)
   if (knownConstitution) then
     write(fileunit,'(a)') '(constitution)'//char(9)//trim(phase_constitution(p))
     do e = 1,phase_Noutput(p)
       write(fileunit,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
     enddo
   endif
 enddo

 close(fileunit)

! allocate memory for state management

 allocate(constitutive_state0(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_partionedState0(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_subState0(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_state(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_dotState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_previousDotState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_previousDotState2(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_relevantState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems))
 allocate(constitutive_sizeDotState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;   constitutive_sizeDotState = 0_pInt
 allocate(constitutive_sizeState(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)) ;      constitutive_sizeState = 0_pInt
 allocate(constitutive_sizePostResults(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems)); constitutive_sizePostResults = 0_pInt
 
 do e = 1,mesh_NcpElems                                  ! loop over elements
   myNgrains = homogenization_Ngrains(mesh_element(3,e)) 
   do i = 1,FE_Nips(mesh_element(2,e))                   ! loop over IPs
     do g = 1,myNgrains                                  ! loop over grains
       selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g)
       myInstance = phase_constitutionInstance(material_phase(g,i,e))
       select case(phase_constitution(material_phase(g,i,e)))  
       
         case (constitutive_j2_label)
           allocate(constitutive_state0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_relevantState(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           constitutive_state0(g,i,e)%p =           constitutive_j2_stateInit(myInstance)
           constitutive_relevantState(g,i,e)%p =    constitutive_j2_relevantState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_j2_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_j2_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_j2_sizePostResults(myInstance)
           
         case (constitutive_phenopowerlaw_label)
           allocate(constitutive_state0(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
           allocate(constitutive_relevantState(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
           constitutive_state0(g,i,e)%p =           constitutive_phenopowerlaw_stateInit(myInstance)
           constitutive_relevantState(g,i,e)%p =    constitutive_phenopowerlaw_relevantState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_phenopowerlaw_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_phenopowerlaw_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_phenopowerlaw_sizePostResults(myInstance)
           
         case (constitutive_dislotwin_label)
           allocate(constitutive_state0(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
           allocate(constitutive_relevantState(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
           constitutive_state0(g,i,e)%p =           constitutive_dislotwin_stateInit(myInstance)
           constitutive_relevantState(g,i,e)%p =    constitutive_dislotwin_relevantState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_dislotwin_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_dislotwin_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_dislotwin_sizePostResults(myInstance)
           
         case (constitutive_nonlocal_label)
           allocate(constitutive_state0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_relevantState(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           constitutive_state0(g,i,e)%p =           constitutive_nonlocal_stateInit(myInstance)
           constitutive_relevantState(g,i,e)%p =    constitutive_nonlocal_relevantState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_nonlocal_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_nonlocal_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_nonlocal_sizePostResults(myInstance)
           
         case default
           call IO_error(200,material_phase(g,i,e))      ! unknown constitution
           
       end select
       constitutive_partionedState0(g,i,e)%p =  constitutive_state0(g,i,e)%p
     enddo
   enddo
 enddo
 
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_maxSizeDotState    = maxval(constitutive_sizeDotState)
 constitutive_maxSizePostResults = maxval(constitutive_sizePostResults)

 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '<<<+-  constitutive init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'constitutive_state0:          ', shape(constitutive_state0)
 write(6,'(a32,x,7(i5,x))') 'constitutive_partionedState0: ', shape(constitutive_partionedState0)
 write(6,'(a32,x,7(i5,x))') 'constitutive_subState0:       ', shape(constitutive_subState0)
 write(6,'(a32,x,7(i5,x))') 'constitutive_state:           ', shape(constitutive_state)
 write(6,'(a32,x,7(i5,x))') 'constitutive_relevantState:   ', shape(constitutive_relevantState)
 write(6,'(a32,x,7(i5,x))') 'constitutive_dotState:        ', shape(constitutive_dotState)
 write(6,'(a32,x,7(i5,x))') 'constitutive_previousDotState:', shape(constitutive_previousDotState)
 write(6,'(a32,x,7(i5,x))') 'constitutive_sizeState:       ', shape(constitutive_sizeState)
 write(6,'(a32,x,7(i5,x))') 'constitutive_sizeDotState:    ', shape(constitutive_sizeDotState)
 write(6,'(a32,x,7(i5,x))') 'constitutive_sizePostResults: ', shape(constitutive_sizePostResults)
 write(6,*)
 write(6,'(a32,x,7(i5,x))') 'maxSizeState:       ', constitutive_maxSizeState
 write(6,'(a32,x,7(i5,x))') 'maxSizeDotState:    ', constitutive_maxSizeDotState
 write(6,'(a32,x,7(i5,x))') 'maxSizePostResults: ', constitutive_maxSizePostResults
 call flush(6)
 !$OMP END CRITICAL (write2out)
 return

endsubroutine


function constitutive_homogenizedC(ipc,ip,el)
!*********************************************************************
!* This function returns the homogenized elacticity matrix           *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal
 implicit none

 !* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal), dimension(6,6) :: constitutive_homogenizedC

 select case (phase_constitution(material_phase(ipc,ip,el)))
 
   case (constitutive_j2_label)
     constitutive_homogenizedC = constitutive_j2_homogenizedC(constitutive_state,ipc,ip,el)
     
   case (constitutive_phenopowerlaw_label)
     constitutive_homogenizedC = constitutive_phenopowerlaw_homogenizedC(constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     constitutive_homogenizedC = constitutive_dislotwin_homogenizedC(constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     constitutive_homogenizedC = constitutive_nonlocal_homogenizedC(constitutive_state,ipc,ip,el)
     
 end select

return
endfunction

function constitutive_averageBurgers(ipc,ip,el)
!*********************************************************************
!* This function returns the average length of Burgers vector        *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal
 implicit none

 !* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal) :: constitutive_averageBurgers

 select case (phase_constitution(material_phase(ipc,ip,el)))
 
   case (constitutive_j2_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_j2_averageBurgers(constitutive_state,ipc,ip,el)
     
   case (constitutive_phenopowerlaw_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_phenopowerlaw_averageBurgers(constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_dislotwin_averageBurgers(constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_nonlocal_averageBurgers(constitutive_state,ipc,ip,el)
     
 end select

return
endfunction


subroutine constitutive_microstructure(Temperature,Tstar_v,Fe,Fp,ipc,ip,el)
!*********************************************************************
!* This function calculates from state needed variables              *
!* INPUT:                                                            *
!*  - state           : state variables                              *
!*  - Tp              : temperature                                  *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec,      only: pReal,pInt
 use material,  only: phase_constitution, &
                      material_phase, &
                      homogenization_maxNgrains
 use mesh,      only: mesh_NcpElems, &
                      mesh_maxNips
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal
 implicit none

!* Definition of variables
integer(pInt), intent(in) :: ipc,ip,el
real(pReal), intent(in) :: Temperature
real(pReal), dimension(6) :: Tstar_v
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: Fe, Fp

 select case (phase_constitution(material_phase(ipc,ip,el)))
 
   case (constitutive_j2_label)
     call constitutive_j2_microstructure(Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_phenopowerlaw_label)
     call constitutive_phenopowerlaw_microstructure(Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     call constitutive_dislotwin_microstructure(Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     call constitutive_nonlocal_microstructure(constitutive_state, Temperature, Tstar_v, Fe, Fp, ipc, ip, el)
     
 end select

endsubroutine


subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, Temperature, ipc, ip, el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - Lp              : plastic velocity gradient                    *
!*  - dLp_dTstar      : derivative of Lp (4th-order tensor)          *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal) Temperature
 real(pReal), dimension(6) :: Tstar_v
 real(pReal), dimension(3,3) :: Lp
 real(pReal), dimension(9,9) :: dLp_dTstar

 select case (phase_constitution(material_phase(ipc,ip,el)))
 
   case (constitutive_j2_label)
     call constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_phenopowerlaw_label)
     call constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     call constitutive_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     call constitutive_nonlocal_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, Temperature, constitutive_state, ipc, ip, el)
     
 end select

 return
endsubroutine


subroutine constitutive_collectDotState(Tstar_v, subTstar0_v, Fe, Fp, Temperature, misorientation, subdt, ipc, ip, el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the rate of change of microstructure                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotState : evolution of state variable            *
!*********************************************************************
use prec, only:     pReal, pInt
use debug, only:    debug_cumDotStateCalls, &
                    debug_cumDotStateTicks
use mesh, only:     mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors
use material, only: phase_constitution, &
                    material_phase, &
                    homogenization_maxNgrains
use constitutive_j2
use constitutive_phenopowerlaw
use constitutive_dislotwin
use constitutive_nonlocal

implicit none

!*** input  variables
integer(pInt), intent(in) ::    ipc, ip, el
real(pReal), intent(in) ::      Temperature, &
                                subdt
real(pReal), dimension(4,mesh_maxNipNeighbors), intent(in) :: &
                                misorientation
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                Fe, &
                                Fp
real(pReal), dimension(6), intent(in) :: &
                                Tstar_v, &
                                subTstar0_v

!*** local variables
integer(pLongInt)               tick, tock, & 
                                tickrate, &
                                maxticks

call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

select case (phase_constitution(material_phase(ipc,ip,el)))

  case (constitutive_j2_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_j2_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
 
  case (constitutive_phenopowerlaw_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_phenopowerlaw_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
 
  case (constitutive_dislotwin_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_dislotwin_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
 
  case (constitutive_nonlocal_label)
    call constitutive_nonlocal_dotState(constitutive_dotState, Tstar_v, subTstar0_v, Fe, Fp, Temperature, misorientation, subdt, &
                                        constitutive_state, constitutive_subState0, subdt, ipc, ip, el)
 
end select

call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
debug_cumDotStateCalls = debug_cumDotStateCalls + 1_pInt
debug_cumDotStateTicks = debug_cumDotStateTicks + tock-tick
if (tock < tick) debug_cumDotStateTicks  = debug_cumDotStateTicks + maxticks

return
endsubroutine


function constitutive_dotTemperature(Tstar_v,Temperature,ipc,ip,el)
!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the rate of change of microstructure                  *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - state           : current microstructure                       *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!* OUTPUT:                                                           *
!*  - constitutive_dotTemperature : evolution of temperature         *
!*********************************************************************
 use prec, only: pReal,pInt
 use material, only: phase_constitution,material_phase
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal
 implicit none

!* Definition of variables
 integer(pInt) ipc,ip,el
 real(pReal) Temperature
 real(pReal) constitutive_dotTemperature
 real(pReal), dimension(6) :: Tstar_v

 select case (phase_constitution(material_phase(ipc,ip,el)))
 
   case (constitutive_j2_label)
     constitutive_dotTemperature = constitutive_j2_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_phenopowerlaw_label)
     constitutive_dotTemperature = constitutive_phenopowerlaw_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     constitutive_dotTemperature = constitutive_dislotwin_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     constitutive_dotTemperature = constitutive_nonlocal_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
     
 end select
 return
endfunction


function constitutive_postResults(Tstar_v, subTstar0_v, Fe, Fp, Temperature, misorientation, dt, subdt, ipc, ip, el)
!*********************************************************************
!* return array of constitutive results                              *
!* INPUT:                                                            *
!*  - Tstar_v         : 2nd Piola Kirchhoff stress tensor (Mandel)   *
!*  - dt              : current time increment                       *
!*  - ipc             : component-ID of current integration point    *
!*  - ip              : current integration point                    *
!*  - el              : current element                              *
!*********************************************************************
 use prec, only:     pReal,pInt
 use mesh, only:     mesh_NcpElems, &
                     mesh_maxNips, &
                     mesh_maxNipNeighbors
 use material, only: phase_constitution, &
                     material_phase, &
                     homogenization_maxNgrains
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_dislotwin
 use constitutive_nonlocal
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: ipc,ip,el
 real(pReal), intent(in) :: dt, Temperature, subdt
 real(pReal), dimension(6), intent(in) :: Tstar_v, subTstar0_v
 real(pReal), dimension(4,mesh_maxNipNeighbors), intent(in) :: misorientation
 real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: Fe, Fp
 real(pReal), dimension(constitutive_sizePostResults(ipc,ip,el)) :: constitutive_postResults

 constitutive_postResults = 0.0_pReal
 select case (phase_constitution(material_phase(ipc,ip,el)))
 
   case (constitutive_j2_label)
     constitutive_postResults = constitutive_j2_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
     
   case (constitutive_phenopowerlaw_label)
     constitutive_postResults = constitutive_phenopowerlaw_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     constitutive_postResults = constitutive_dislotwin_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     constitutive_postResults = constitutive_nonlocal_postResults(Tstar_v, subTstar0_v, Fe, Fp, Temperature, misorientation, &
                                                                  dt, subdt, constitutive_state, constitutive_subState0, &
                                                                  constitutive_dotstate, ipc, ip, el)
 end select
 
return

endfunction

END MODULE