! Copyright 2011 Max-Planck-Institut für Eisenforschung GmbH
!
! This file is part of DAMASK,
! the Düsseldorf Advanced MAterial Simulation Kit.
!
! DAMASK is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! DAMASK is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with DAMASK. If not, see <http://www.gnu.org/licenses/>.
!
!##############################################################
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
                                               constitutive_state_backup, &    ! pointer array to backed up microstructure (end of converged time step)
                                               constitutive_dotState, &        ! pointer array to evolution of current microstructure
                                               constitutive_previousDotState,& ! pointer array to previous evolution of current microstructure
                                               constitutive_previousDotState2,&! pointer array to 2nd previous evolution of current microstructure
                                               constitutive_dotState_backup, & ! pointer array to backed up evolution of current microstructure
                                               constitutive_RK4dotState, &     ! pointer array to evolution of microstructure defined by classical Runge-Kutta method
                                               constitutive_aTolState          ! pointer array to absolute state tolerance
type(p_vec), dimension(:,:,:,:), allocatable :: constitutive_RKCK45dotState     ! pointer array to evolution of microstructure used by Cash-Karp Runge-Kutta method
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


!**************************************
!*      Module initialization         *
!**************************************
subroutine constitutive_init()
use prec, only: pReal,pInt
use debug, only: debug_verbosity, debug_selectiveDebugger, debug_e, debug_i, debug_g
use numerics, only: numerics_integrator
use IO, only: IO_error, IO_open_file, IO_open_jobFile
use mesh, only: mesh_maxNips,mesh_NcpElems,mesh_element,FE_Nips
use material
use constitutive_j2
use constitutive_phenopowerlaw
use constitutive_titanmod
use constitutive_dislotwin
use constitutive_nonlocal

implicit none

integer(pInt), parameter :: fileunit = 200
integer(pInt)   g, &                          ! grain number
                i, &                          ! integration point number
                e, &                          ! element number
                gMax, &                       ! maximum number of grains
                iMax, &                       ! maximum number of integration points
                eMax, &                       ! maximum number of elements
                p, &
                s, &
                myInstance,& 
                myNgrains
integer(pInt), dimension(:,:), pointer :: thisSize
character(len=64), dimension(:,:), pointer :: thisOutput
logical knownConstitution


! --- PARSE CONSTITUTIONS FROM CONFIG FILE ---

if(.not. IO_open_file(fileunit,material_configFile)) call IO_error (100) ! corrupt config file
call constitutive_j2_init(fileunit)
call constitutive_phenopowerlaw_init(fileunit)
call constitutive_titanmod_init(fileunit)
call constitutive_dislotwin_init(fileunit)
call constitutive_nonlocal_init(fileunit)  
close(fileunit)


! --- WRITE DESCRIPTION FILE FOR CONSTITUTIVE PHASE OUTPUT ---

if(.not. IO_open_jobFile(fileunit,'outputConstitutive')) then ! problems in writing file
  call IO_error (50) 
endif
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
    case (constitutive_titanmod_label)
      thisOutput => constitutive_titanmod_output
      thisSize   => constitutive_titanmod_sizePostResult
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


! --- ALLOCATION OF STATES ---

gMax = homogenization_maxNgrains
iMax = mesh_maxNips
eMax = mesh_NcpElems

allocate(constitutive_state0(gMax,iMax,eMax))
allocate(constitutive_partionedState0(gMax,iMax,eMax))
allocate(constitutive_subState0(gMax,iMax,eMax))
allocate(constitutive_state(gMax,iMax,eMax))
allocate(constitutive_state_backup(gMax,iMax,eMax))
allocate(constitutive_dotState(gMax,iMax,eMax))
allocate(constitutive_dotState_backup(gMax,iMax,eMax))
allocate(constitutive_aTolState(gMax,iMax,eMax))
allocate(constitutive_sizeDotState(gMax,iMax,eMax)) ;          constitutive_sizeDotState = 0_pInt
allocate(constitutive_sizeState(gMax,iMax,eMax)) ;                constitutive_sizeState = 0_pInt
allocate(constitutive_sizePostResults(gMax,iMax,eMax));     constitutive_sizePostResults = 0_pInt
if (any(numerics_integrator == 1)) then
  allocate(constitutive_previousDotState(gMax,iMax,eMax))
  allocate(constitutive_previousDotState2(gMax,iMax,eMax))
endif
if (any(numerics_integrator == 4)) then
  allocate(constitutive_RK4dotState(gMax,iMax,eMax)) 
endif
if (any(numerics_integrator == 5)) then
  allocate(constitutive_RKCK45dotState(6,gMax,iMax,eMax))
endif

!$OMP PARALLEL DO PRIVATE(myNgrains,myInstance)
  do e = 1,mesh_NcpElems                                  ! loop over elements
    myNgrains = homogenization_Ngrains(mesh_element(3,e)) 
    do i = 1,FE_Nips(mesh_element(2,e))                   ! loop over IPs
      do g = 1,myNgrains                                  ! loop over grains
        myInstance = phase_constitutionInstance(material_phase(g,i,e))
        select case(phase_constitution(material_phase(g,i,e)))  
        
          case (constitutive_j2_label)
            allocate(constitutive_state0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
            allocate(constitutive_partionedState0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
            allocate(constitutive_subState0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
            allocate(constitutive_state(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
            allocate(constitutive_state_backup(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
            allocate(constitutive_aTolState(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
            allocate(constitutive_dotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
            allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
            if (any(numerics_integrator == 1)) then
              allocate(constitutive_previousDotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
              allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
            endif
            if (any(numerics_integrator == 4)) then
              allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance))) 
            endif
            if (any(numerics_integrator == 5)) then
              do s = 1,6
                allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_j2_sizeDotState(myInstance))) 
              enddo
            endif
            constitutive_state0(g,i,e)%p =           constitutive_j2_stateInit(myInstance)
            constitutive_aTolState(g,i,e)%p =        constitutive_j2_aTolState(myInstance)
            constitutive_sizeState(g,i,e) =          constitutive_j2_sizeState(myInstance)
            constitutive_sizeDotState(g,i,e) =       constitutive_j2_sizeDotState(myInstance)
            constitutive_sizePostResults(g,i,e) =    constitutive_j2_sizePostResults(myInstance)
           
          case (constitutive_phenopowerlaw_label)
            allocate(constitutive_state0(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
            allocate(constitutive_partionedState0(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
            allocate(constitutive_subState0(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
            allocate(constitutive_state(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
            allocate(constitutive_state_backup(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
            allocate(constitutive_aTolState(g,i,e)%p(constitutive_phenopowerlaw_sizeState(myInstance)))
            allocate(constitutive_dotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
            allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
            if (any(numerics_integrator == 1)) then
              allocate(constitutive_previousDotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
              allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
            endif
            if (any(numerics_integrator == 4)) then
              allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance))) 
            endif
            if (any(numerics_integrator == 5)) then
              do s = 1,6
                allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance))) 
              enddo
            endif
            constitutive_state0(g,i,e)%p =           constitutive_phenopowerlaw_stateInit(myInstance)
            constitutive_aTolState(g,i,e)%p =        constitutive_phenopowerlaw_aTolState(myInstance)
            constitutive_sizeState(g,i,e) =          constitutive_phenopowerlaw_sizeState(myInstance)
            constitutive_sizeDotState(g,i,e) =       constitutive_phenopowerlaw_sizeDotState(myInstance)
            constitutive_sizePostResults(g,i,e) =    constitutive_phenopowerlaw_sizePostResults(myInstance)
            
          case (constitutive_titanmod_label)
            allocate(constitutive_state0(g,i,e)%p(constitutive_titanmod_sizeState(myInstance)))
            allocate(constitutive_partionedState0(g,i,e)%p(constitutive_titanmod_sizeState(myInstance)))
            allocate(constitutive_subState0(g,i,e)%p(constitutive_titanmod_sizeState(myInstance)))
            allocate(constitutive_state(g,i,e)%p(constitutive_titanmod_sizeState(myInstance)))
            allocate(constitutive_state_backup(g,i,e)%p(constitutive_titanmod_sizeState(myInstance)))
            allocate(constitutive_aTolState(g,i,e)%p(constitutive_titanmod_sizeState(myInstance)))
            allocate(constitutive_dotState(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
            allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
            if (any(numerics_integrator == 1)) then
              allocate(constitutive_previousDotState(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
              allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
            endif
            if (any(numerics_integrator == 4)) then
              allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance))) 
            endif
            if (any(numerics_integrator == 5)) then
              do s = 1,6
                allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance))) 
              enddo
            endif
            constitutive_state0(g,i,e)%p =           constitutive_titanmod_stateInit(myInstance)
            constitutive_aTolState(g,i,e)%p =        constitutive_titanmod_aTolState(myInstance)
            constitutive_sizeState(g,i,e) =          constitutive_titanmod_sizeState(myInstance)
            constitutive_sizeDotState(g,i,e) =       constitutive_titanmod_sizeDotState(myInstance)
            constitutive_sizePostResults(g,i,e) =    constitutive_titanmod_sizePostResults(myInstance)
          
          case (constitutive_dislotwin_label)
            allocate(constitutive_state0(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
            allocate(constitutive_partionedState0(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
            allocate(constitutive_subState0(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
            allocate(constitutive_state(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
            allocate(constitutive_state_backup(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
            allocate(constitutive_aTolState(g,i,e)%p(constitutive_dislotwin_sizeState(myInstance)))
            allocate(constitutive_dotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
            allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
            if (any(numerics_integrator == 1)) then
              allocate(constitutive_previousDotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
              allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
            endif
            if (any(numerics_integrator == 4)) then
              allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance))) 
            endif
            if (any(numerics_integrator == 5)) then
              do s = 1,6
                allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance))) 
              enddo
            endif
            constitutive_state0(g,i,e)%p =           constitutive_dislotwin_stateInit(myInstance)
            constitutive_aTolState(g,i,e)%p =        constitutive_dislotwin_aTolState(myInstance)
            constitutive_sizeState(g,i,e) =          constitutive_dislotwin_sizeState(myInstance)
            constitutive_sizeDotState(g,i,e) =       constitutive_dislotwin_sizeDotState(myInstance)
            constitutive_sizePostResults(g,i,e) =    constitutive_dislotwin_sizePostResults(myInstance)
            
          case (constitutive_nonlocal_label)
            allocate(constitutive_state0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
            allocate(constitutive_partionedState0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
            allocate(constitutive_subState0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
            allocate(constitutive_state(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
            allocate(constitutive_state_backup(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
            allocate(constitutive_aTolState(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
            allocate(constitutive_dotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
            allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
            if (any(numerics_integrator == 1)) then
              allocate(constitutive_previousDotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
              allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
            endif
            if (any(numerics_integrator == 4)) then
              allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance))) 
            endif
            if (any(numerics_integrator == 5)) then
              do s = 1,6
                allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance))) 
              enddo
            endif
            constitutive_state0(g,i,e)%p =           constitutive_nonlocal_stateInit(myInstance)
            constitutive_aTolState(g,i,e)%p =        constitutive_nonlocal_aTolState(myInstance)
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
!$OMP END PARALLEL DO

constitutive_maxSizeState       = maxval(constitutive_sizeState)
constitutive_maxSizeDotState    = maxval(constitutive_sizeDotState)
constitutive_maxSizePostResults = maxval(constitutive_sizePostResults)

!$OMP CRITICAL (write2out)
  write(6,*)
  write(6,*) '<<<+-  constitutive init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
  if (debug_verbosity > 0) then
    write(6,'(a32,x,7(i8,x))') 'constitutive_state0:          ', shape(constitutive_state0)
    write(6,'(a32,x,7(i8,x))') 'constitutive_partionedState0: ', shape(constitutive_partionedState0)
    write(6,'(a32,x,7(i8,x))') 'constitutive_subState0:       ', shape(constitutive_subState0)
    write(6,'(a32,x,7(i8,x))') 'constitutive_state:           ', shape(constitutive_state)
    write(6,'(a32,x,7(i8,x))') 'constitutive_aTolState:       ', shape(constitutive_aTolState)
    write(6,'(a32,x,7(i8,x))') 'constitutive_dotState:        ', shape(constitutive_dotState)
    write(6,'(a32,x,7(i8,x))') 'constitutive_sizeState:       ', shape(constitutive_sizeState)
    write(6,'(a32,x,7(i8,x))') 'constitutive_sizeDotState:    ', shape(constitutive_sizeDotState)
    write(6,'(a32,x,7(i8,x))') 'constitutive_sizePostResults: ', shape(constitutive_sizePostResults)
    write(6,*)
    write(6,'(a32,x,7(i8,x))') 'maxSizeState:       ', constitutive_maxSizeState
    write(6,'(a32,x,7(i8,x))') 'maxSizeDotState:    ', constitutive_maxSizeDotState
    write(6,'(a32,x,7(i8,x))') 'maxSizePostResults: ', constitutive_maxSizePostResults
  endif
  call flush(6)
!$OMP END CRITICAL (write2out)

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
 use constitutive_titanmod
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

   case (constitutive_titanmod_label)
     constitutive_homogenizedC = constitutive_titanmod_homogenizedC(constitutive_state,ipc,ip,el)
 
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
 use constitutive_titanmod
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

   case (constitutive_titanmod_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_titanmod_averageBurgers(constitutive_state,ipc,ip,el)
     
   case (constitutive_dislotwin_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_dislotwin_averageBurgers(constitutive_state,ipc,ip,el)
     
   case (constitutive_nonlocal_label)
     constitutive_averageBurgers = 2.5e-10_pReal !constitutive_nonlocal_averageBurgers(constitutive_state,ipc,ip,el)
     
 end select

return
endfunction



!*********************************************************************
!* This function calculates from state needed variables              *
!*********************************************************************
subroutine constitutive_microstructure(Temperature, Fe, ipc, ip, el)
use prec,      only: pReal,pInt
use material,  only: phase_constitution, &
                     material_phase, &
                     homogenization_maxNgrains
use mesh,      only: mesh_NcpElems, &
                     mesh_maxNips, &
                     mesh_maxNipNeighbors
use constitutive_j2,            only: constitutive_j2_label, &
                                      constitutive_j2_microstructure
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_label, &
                                      constitutive_phenopowerlaw_microstructure
use constitutive_titanmod,      only: constitutive_titanmod_label, &
                                      constitutive_titanmod_microstructure
use constitutive_dislotwin,     only: constitutive_dislotwin_label, &
                                      constitutive_dislotwin_microstructure
use constitutive_nonlocal,      only: constitutive_nonlocal_label, &
                                      constitutive_nonlocal_microstructure
implicit none

!*** input variables ***!
integer(pInt), intent(in)::               ipc, &        ! component-ID of current integration point
                                          ip, &         ! current integration point
                                          el            ! current element
real(pReal), intent(in) ::                Temperature
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                          Fe            ! elastic deformation gradient

!*** output variables ***!

!*** local variables ***!


select case (phase_constitution(material_phase(ipc,ip,el)))
 
  case (constitutive_j2_label)
    call constitutive_j2_microstructure(Temperature,constitutive_state,ipc,ip,el)
     
  case (constitutive_phenopowerlaw_label)
    call constitutive_phenopowerlaw_microstructure(Temperature,constitutive_state,ipc,ip,el)
  
  case (constitutive_titanmod_label)
    call constitutive_titanmod_microstructure(Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_dislotwin_label)
    call constitutive_dislotwin_microstructure(Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_nonlocal_label)
    call constitutive_nonlocal_microstructure(constitutive_state, Temperature, Fe, ipc, ip, el)
     
end select

endsubroutine



!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the velocity gradient                                 *
!*********************************************************************
subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, Temperature, ipc, ip, el)

use prec, only: pReal,pInt
use material, only: phase_constitution, &
                    material_phase
use constitutive_j2,            only: constitutive_j2_label, &
                                      constitutive_j2_LpAndItsTangent
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_label, &
                                      constitutive_phenopowerlaw_LpAndItsTangent
use constitutive_titanmod,      only: constitutive_titanmod_label, &
                                      constitutive_titanmod_LpAndItsTangent
use constitutive_dislotwin,     only: constitutive_dislotwin_label, &
                                      constitutive_dislotwin_LpAndItsTangent
use constitutive_nonlocal,      only: constitutive_nonlocal_label, &
                                      constitutive_nonlocal_LpAndItsTangent
implicit none


!*** input variables ***!
integer(pInt), intent(in)::                 ipc, &        ! component-ID of current integration point
                                            ip, &         ! current integration point
                                            el            ! current element
real(pReal), intent(in) ::                  Temperature
real(pReal), dimension(6), intent(in) ::    Tstar_v       ! 2nd Piola-Kirchhoff stress

!*** output variables ***!
real(pReal), dimension(3,3), intent(out) :: Lp            ! plastic velocity gradient
real(pReal), dimension(9,9), intent(out) :: dLp_dTstar    ! derivative of Lp with respect to Tstar (4th-order tensor)


!*** local variables ***!


select case (phase_constitution(material_phase(ipc,ip,el)))

  case (constitutive_j2_label)
    call constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_phenopowerlaw_label)
    call constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
  
  case (constitutive_titanmod_label)
    call constitutive_titanmod_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_dislotwin_label)
    call constitutive_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_nonlocal_label)
    call constitutive_nonlocal_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, Temperature, constitutive_state, ipc, ip, el)
   
end select

endsubroutine



!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the rate of change of microstructure                  *
!*********************************************************************
subroutine constitutive_collectDotState(Tstar_v, Fe, Fp, Temperature, subdt, orientation, ipc, ip, el)

use prec, only:     pReal, pInt
use debug, only:    debug_cumDotStateCalls, &
                    debug_cumDotStateTicks, &
                    debug_verbosity
use mesh, only:     mesh_NcpElems, &
                    mesh_maxNips, &
                    mesh_maxNipNeighbors
use material, only: phase_constitution, &
                    material_phase, &
                    homogenization_maxNgrains
use constitutive_j2, only:            constitutive_j2_dotState, &
                                      constitutive_j2_label
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_dotState, &
                                      constitutive_phenopowerlaw_label
use constitutive_titanmod, only:      constitutive_titanmod_dotState, &
                                      constitutive_titanmod_label
use constitutive_dislotwin, only:     constitutive_dislotwin_dotState, &
                                      constitutive_dislotwin_label
use constitutive_nonlocal, only:      constitutive_nonlocal_dotState, &
                                      constitutive_nonlocal_label

implicit none

!*** input  variables
integer(pInt), intent(in) ::    ipc, &        ! component-ID of current integration point
                                ip, &         ! current integration point
                                el            ! current element
real(pReal), intent(in) ::      Temperature, &
                                subdt         ! timestep
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                Fe, &         ! elastic deformation gradient
                                Fp            ! plastic deformation gradient
real(pReal), dimension(4,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems), intent(in) :: &
                                orientation   ! crystal orientation (quaternion)
real(pReal), dimension(6), intent(in) :: &
                                Tstar_v       ! 2nd Piola Kirchhoff stress tensor (Mandel)

!*** output variables ***!

!*** local variables
integer(pLongInt)               tick, tock, & 
                                tickrate, &
                                maxticks

if (debug_verbosity > 0) then
  call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
endif

select case (phase_constitution(material_phase(ipc,ip,el)))

  case (constitutive_j2_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_j2_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
 
  case (constitutive_phenopowerlaw_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_phenopowerlaw_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)

  case (constitutive_titanmod_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_titanmod_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
  
  case (constitutive_dislotwin_label)
    constitutive_dotState(ipc,ip,el)%p = constitutive_dislotwin_dotState(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
 
  case (constitutive_nonlocal_label)
    call constitutive_nonlocal_dotState(constitutive_dotState(ipc,ip,el), Tstar_v, Fe, Fp, Temperature, constitutive_state, &
                                        constitutive_aTolState, subdt, orientation, ipc, ip, el)
 
end select

if (debug_verbosity > 0) then
  call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
  !$OMP CRITICAL (debugTimingDotState)
    debug_cumDotStateCalls = debug_cumDotStateCalls + 1_pInt
    debug_cumDotStateTicks = debug_cumDotStateTicks + tock-tick
    !$OMP FLUSH (debug_cumDotStateTicks)
    if (tock < tick) debug_cumDotStateTicks  = debug_cumDotStateTicks + maxticks
  !$OMP END CRITICAL (debugTimingDotState)
endif

endsubroutine



!*********************************************************************
!* This subroutine contains the constitutive equation for            *
!* calculating the rate of change of microstructure                  *
!*********************************************************************
function constitutive_dotTemperature(Tstar_v,Temperature,ipc,ip,el)

use prec, only:     pReal,pInt
use debug, only:    debug_cumDotTemperatureCalls, &
                    debug_cumDotTemperatureTicks, &
                    debug_verbosity
use material, only: phase_constitution, &
                    material_phase
use constitutive_j2, only:            constitutive_j2_dotTemperature, &
                                      constitutive_j2_label
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_dotTemperature, &
                                      constitutive_phenopowerlaw_label
use constitutive_titanmod, only:      constitutive_titanmod_dotTemperature, &
                                      constitutive_titanmod_label
use constitutive_dislotwin, only:     constitutive_dislotwin_dotTemperature, &
                                      constitutive_dislotwin_label
use constitutive_nonlocal, only:      constitutive_nonlocal_dotTemperature, &
                                      constitutive_nonlocal_label
implicit none

!*** input  variables
integer(pInt), intent(in) ::    ipc, &        ! component-ID of current integration point
                                ip, &         ! current integration point
                                el            ! current element
real(pReal), intent(in) ::      Temperature
real(pReal), dimension(6), intent(in) :: &
                                Tstar_v       ! 2nd Piola Kirchhoff stress tensor (Mandel)

!*** output variables ***!
real(pReal)                     constitutive_dotTemperature   ! evolution of temperature

!*** local variables
integer(pLongInt)               tick, tock, & 
                                tickrate, &
                                maxticks


if (debug_verbosity > 0) then
  call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
endif

select case (phase_constitution(material_phase(ipc,ip,el)))

  case (constitutive_j2_label)
    constitutive_dotTemperature = constitutive_j2_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_phenopowerlaw_label)
    constitutive_dotTemperature = constitutive_phenopowerlaw_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)

  case (constitutive_titanmod_label)
    constitutive_dotTemperature = constitutive_titanmod_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_dislotwin_label)
    constitutive_dotTemperature = constitutive_dislotwin_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
  case (constitutive_nonlocal_label)
    constitutive_dotTemperature = constitutive_nonlocal_dotTemperature(Tstar_v,Temperature,constitutive_state,ipc,ip,el)
   
end select

if (debug_verbosity > 0) then
  call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
  !$OMP CRITICAL (debugTimingDotTemperature)
    debug_cumDotTemperatureCalls = debug_cumDotTemperatureCalls + 1_pInt
    debug_cumDotTemperatureTicks = debug_cumDotTemperatureTicks + tock-tick
    !$OMP FLUSH (debug_cumDotTemperatureTicks)
    if (tock < tick) debug_cumDotTemperatureTicks  = debug_cumDotTemperatureTicks + maxticks
  !$OMP END CRITICAL (debugTimingDotTemperature)
endif

endfunction



function constitutive_postResults(Tstar_v, Fe, Temperature, dt, ipc, ip, el)
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
use constitutive_j2, only:            constitutive_j2_postResults, &
                                      constitutive_j2_label
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_postResults, &
                                      constitutive_phenopowerlaw_label
use constitutive_titanmod, only:      constitutive_titanmod_postResults, &
                                      constitutive_titanmod_label
use constitutive_dislotwin, only:     constitutive_dislotwin_postResults, &
                                      constitutive_dislotwin_label
use constitutive_nonlocal, only:      constitutive_nonlocal_postResults, &
                                      constitutive_nonlocal_label
implicit none

!*** input  variables
integer(pInt), intent(in) ::    ipc, &        ! component-ID of current integration point
                                ip, &         ! current integration point
                                el            ! current element
real(pReal), intent(in) ::      Temperature, &
                                dt            ! timestep
real(pReal), dimension(3,3), intent(in) :: &
                                Fe            ! elastic deformation gradient
real(pReal), dimension(6), intent(in) :: &
                                Tstar_v       ! 2nd Piola Kirchhoff stress tensor (Mandel)

!*** output variables ***!
real(pReal), dimension(constitutive_sizePostResults(ipc,ip,el)) :: constitutive_postResults

!*** local variables


constitutive_postResults = 0.0_pReal
select case (phase_constitution(material_phase(ipc,ip,el)))

  case (constitutive_j2_label)
    constitutive_postResults = constitutive_j2_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
   
  case (constitutive_phenopowerlaw_label)
    constitutive_postResults = constitutive_phenopowerlaw_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
  
  case (constitutive_titanmod_label)
    constitutive_postResults = constitutive_titanmod_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
   
  case (constitutive_dislotwin_label)
    constitutive_postResults = constitutive_dislotwin_postResults(Tstar_v,Temperature,dt,constitutive_state,ipc,ip,el)
   
  case (constitutive_nonlocal_label)
    constitutive_postResults = constitutive_nonlocal_postResults(Tstar_v, Fe, Temperature, dt, constitutive_state, &
                                                                 constitutive_dotstate, ipc, ip, el)
end select
 
endfunction

END MODULE