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
!--------------------------------------------------------------------------------------------------
!* $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!! Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief elasticity, plasticity, internal microstructure state
!--------------------------------------------------------------------------------------------------
module constitutive
 use prec, only: &
   pInt, &
   pReal, &
   p_vec
 
 implicit none
 private
 type(p_vec),   public, dimension(:,:,:), allocatable :: &
   constitutive_state0, &                                                                            !< pointer array to microstructure at start of FE inc
   constitutive_partionedState0, &                                                                   !< pointer array to microstructure at start of homogenization inc
   constitutive_subState0, &                                                                         !< pointer array to microstructure at start of crystallite inc
   constitutive_state, &                                                                             !< pointer array to current microstructure (end of converged time step)
   constitutive_state_backup, &                                                                      !< pointer array to backed up microstructure (end of converged time step)
   constitutive_dotState, &                                                                          !< pointer array to evolution of current microstructure
   constitutive_deltaState, &                                                                        !< pointer array to incremental change of current microstructure
   constitutive_previousDotState,&                                                                   !< pointer array to previous evolution of current microstructure
   constitutive_previousDotState2,&                                                                  !< pointer array to 2nd previous evolution of current microstructure
   constitutive_dotState_backup, &                                                                   !< pointer array to backed up evolution of current microstructure
   constitutive_RK4dotState, &                                                                       !< pointer array to evolution of microstructure defined by classical Runge-Kutta method
   constitutive_aTolState                                                                            !< pointer array to absolute state tolerance
 type(p_vec),   public, dimension(:,:,:,:), allocatable :: &
   constitutive_RKCK45dotState                                                                       !< pointer array to evolution of microstructure used by Cash-Karp Runge-Kutta method
 integer(pInt), public, dimension(:,:,:), allocatable :: &
   constitutive_sizeDotState, &                                                                      !< size of dotState array
   constitutive_sizeState, &                                                                         !< size of state array per grain
   constitutive_sizePostResults                                                                      !< size of postResults array per grain
 integer(pInt), public :: &
   constitutive_maxSizeDotState, &
   constitutive_maxSizePostResults
 integer(pInt), private :: &
   constitutive_maxSizeState
  
 character (len=*), parameter, public :: constitutive_hooke_label = 'hooke'
 
 public :: & 
   constitutive_init, &
   constitutive_homogenizedC, &
   constitutive_averageBurgers, &
   constitutive_microstructure, &
   constitutive_LpAndItsTangent, &
   constitutive_TandItsTangent, &
   constitutive_collectDotState, &
   constitutive_collectDeltaState, &
   constitutive_dotTemperature, &
   constitutive_postResults
 
 private :: &
   constitutive_hooke_TandItsTangent
 
contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates arrays pointing to array of the various constitutive modules
!--------------------------------------------------------------------------------------------------
subroutine constitutive_init
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only: &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use numerics, only: &
   numerics_integrator
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_write_jobBinaryIntFile
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, &
   mesh_ipNeighborhood, &
   mesh_maxNipNeighbors, &
   FE_Nips, &
   FE_NipNeighbors, &
   FE_geomtype
 use material, only: &
   material_phase, &
   material_Nphase, &
   material_localFileExt, &    
   material_configFile, &    
   phase_name, &
   phase_elasticity, &
   phase_plasticity, &
   phase_plasticityInstance, &
   phase_Noutput, &
   homogenization_Ngrains, &
   homogenization_maxNgrains
 use constitutive_none
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_titanmod
 use constitutive_dislotwin
 use constitutive_nonlocal

 implicit none
 integer(pInt), parameter :: fileunit = 200_pInt
 integer(pInt) :: &
  g, &                                                                                              !< grain number
  i, &                                                                                              !< integration point number
  e, &                                                                                              !< element number
  gMax, &                                                                                           !< maximum number of grains
  iMax, &                                                                                           !< maximum number of integration points
  eMax, &                                                                                           !< maximum number of elements
  p, &
  s, &
  myInstance,& 
  myNgrains
 integer(pInt), dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 logical :: knownPlasticity
 
 
!--------------------------------------------------------------------------------------------------
! parse plasticities from config file
 if (.not. IO_open_jobFile_stat(fileunit,material_localFileExt)) then                               ! no local material configuration present...
   call IO_open_file(fileunit,material_configFile)                                                  ! ... open material.config file
 endif
 call constitutive_none_init(fileunit)
 call constitutive_j2_init(fileunit)
 call constitutive_phenopowerlaw_init(fileunit)
 call constitutive_titanmod_init(fileunit)
 call constitutive_dislotwin_init(fileunit)
 call constitutive_nonlocal_init(fileunit)  
 close(fileunit)
 
 write(6,*)
 write(6,*) '<<<+-  constitutive init  -+>>>'
 write(6,*) '$Id$'
#include "compilation_info.f90"
 
!--------------------------------------------------------------------------------------------------
! write description file for constitutive phase output
 call IO_write_jobFile(fileunit,'outputConstitutive') 
 do p = 1_pInt,material_Nphase
   i = phase_plasticityInstance(p)                                                                  ! which instance of a plasticity is present phase
   knownPlasticity = .true.                                                                         ! assume valid
   select case(phase_plasticity(p))                                                                 ! split per constitiution
     case (constitutive_none_label)
       thisOutput => NULL()                                                                         ! constitutive_none_output
       thisSize   => NULL()                                                                         ! constitutive_none_sizePostResult
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
       knownPlasticity = .false.
   end select   
   write(fileunit,*)
   write(fileunit,'(a)') '['//trim(phase_name(p))//']'
   write(fileunit,*)
   if (knownPlasticity) then
     write(fileunit,'(a)') '(plasticity)'//char(9)//trim(phase_plasticity(p))
     do e = 1_pInt,phase_Noutput(p)
       write(fileunit,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
     enddo
   endif
 enddo
 close(fileunit)
 
!--------------------------------------------------------------------------------------------------
! allocation of states
 gMax = homogenization_maxNgrains
 iMax = mesh_maxNips
 eMax = mesh_NcpElems
 
 allocate(constitutive_state0(gMax,iMax,eMax))
 allocate(constitutive_partionedState0(gMax,iMax,eMax))
 allocate(constitutive_subState0(gMax,iMax,eMax))
 allocate(constitutive_state(gMax,iMax,eMax))
 allocate(constitutive_state_backup(gMax,iMax,eMax))
 allocate(constitutive_dotState(gMax,iMax,eMax))
 allocate(constitutive_deltaState(gMax,iMax,eMax))
 allocate(constitutive_dotState_backup(gMax,iMax,eMax))
 allocate(constitutive_aTolState(gMax,iMax,eMax))
 allocate(constitutive_sizeDotState(gMax,iMax,eMax)) ;          constitutive_sizeDotState = 0_pInt
 allocate(constitutive_sizeState(gMax,iMax,eMax)) ;                constitutive_sizeState = 0_pInt
 allocate(constitutive_sizePostResults(gMax,iMax,eMax));     constitutive_sizePostResults = 0_pInt
 if (any(numerics_integrator == 1_pInt)) then
   allocate(constitutive_previousDotState(gMax,iMax,eMax))
   allocate(constitutive_previousDotState2(gMax,iMax,eMax))
 endif
 if (any(numerics_integrator == 4_pInt)) then
   allocate(constitutive_RK4dotState(gMax,iMax,eMax)) 
 endif
 if (any(numerics_integrator == 5_pInt)) then
   allocate(constitutive_RKCK45dotState(6,gMax,iMax,eMax))
 endif
 
 do e = 1_pInt,mesh_NcpElems                                                                        ! loop over elements
   myNgrains = homogenization_Ngrains(mesh_element(3,e)) 
   do i = 1_pInt,FE_Nips(FE_geomtype(mesh_element(2,e)))                                            ! loop over IPs
     do g = 1_pInt,myNgrains                                                                        ! loop over grains
       select case(phase_elasticity(material_phase(g,i,e)))                                        
 
         case (constitutive_hooke_label)
                                                                                                    ! valid elasticity but nothing to do             
         case default
           call IO_error(200_pInt,ext_msg=trim(phase_elasticity(material_phase(g,i,e))))            ! unknown elasticity
          
       end select
       myInstance = phase_plasticityInstance(material_phase(g,i,e))
       select case(phase_plasticity(material_phase(g,i,e)))  
       
         case (constitutive_none_label)
           allocate(constitutive_state0(g,i,e)%p(constitutive_none_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_none_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_none_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_none_sizeState(myInstance)))
           allocate(constitutive_state_backup(g,i,e)%p(constitutive_none_sizeState(myInstance)))
           allocate(constitutive_aTolState(g,i,e)%p(constitutive_none_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_none_sizeDotState(myInstance)))
           allocate(constitutive_deltaState(g,i,e)%p(constitutive_none_sizeDotState(myInstance)))
           allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_none_sizeDotState(myInstance)))
           if (any(numerics_integrator == 1_pInt)) then
             allocate(constitutive_previousDotState(g,i,e)%p(constitutive_none_sizeDotState(myInstance)))
             allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_none_sizeDotState(myInstance)))
           endif
           if (any(numerics_integrator == 4_pInt)) then
             allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_none_sizeDotState(myInstance))) 
           endif
           if (any(numerics_integrator == 5_pInt)) then
             do s = 1_pInt,6_pInt
               allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_none_sizeDotState(myInstance))) 
             enddo
           endif
           constitutive_state0(g,i,e)%p =           constitutive_none_stateInit(myInstance)
           constitutive_aTolState(g,i,e)%p =        constitutive_none_aTolState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_none_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_none_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_none_sizePostResults(myInstance)
          
         case (constitutive_j2_label)
           allocate(constitutive_state0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_state_backup(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_aTolState(g,i,e)%p(constitutive_j2_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           allocate(constitutive_deltaState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           if (any(numerics_integrator == 1_pInt)) then
             allocate(constitutive_previousDotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
             allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_j2_sizeDotState(myInstance)))
           endif
           if (any(numerics_integrator == 4_pInt)) then
             allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_j2_sizeDotState(myInstance))) 
           endif
           if (any(numerics_integrator == 5_pInt)) then
             do s = 1_pInt,6_pInt
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
           allocate(constitutive_deltaState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
           allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
           if (any(numerics_integrator == 1_pInt)) then
             allocate(constitutive_previousDotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
             allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance)))
           endif
           if (any(numerics_integrator == 4_pInt)) then
             allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_phenopowerlaw_sizeDotState(myInstance))) 
           endif
           if (any(numerics_integrator == 5_pInt)) then
             do s = 1_pInt,6_pInt
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
           allocate(constitutive_deltaState(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
           allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
           if (any(numerics_integrator == 1_pInt)) then
             allocate(constitutive_previousDotState(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
             allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance)))
           endif
           if (any(numerics_integrator == 4_pInt)) then
             allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_titanmod_sizeDotState(myInstance))) 
           endif
           if (any(numerics_integrator == 5_pInt)) then
             do s = 1_pInt,6_pInt
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
           allocate(constitutive_deltaState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
           allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
           if (any(numerics_integrator == 1_pInt)) then
             allocate(constitutive_previousDotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
             allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance)))
           endif
           if (any(numerics_integrator == 4_pInt)) then
             allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance))) 
           endif
           if (any(numerics_integrator == 5_pInt)) then
             do s = 1_pInt,6_pInt
               allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_dislotwin_sizeDotState(myInstance))) 
             enddo
           endif
           constitutive_state0(g,i,e)%p =           constitutive_dislotwin_stateInit(myInstance)
           constitutive_aTolState(g,i,e)%p =        constitutive_dislotwin_aTolState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_dislotwin_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_dislotwin_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_dislotwin_sizePostResults(myInstance)
           
         case (constitutive_nonlocal_label)
           if(myNgrains/=1_pInt) call IO_error(252_pInt, e,i,g)
           select case(FE_geomtype(mesh_element(2,e)))
             case (7_pInt,8_pInt,9_pInt,10_pInt)
               ! all fine
             case default
               call IO_error(253_pInt,e,i,g)
           end select
           allocate(constitutive_state0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_partionedState0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_subState0(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_state(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_state_backup(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_aTolState(g,i,e)%p(constitutive_nonlocal_sizeState(myInstance)))
           allocate(constitutive_dotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           allocate(constitutive_deltaState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           allocate(constitutive_dotState_backup(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           if (any(numerics_integrator == 1_pInt)) then
             allocate(constitutive_previousDotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
             allocate(constitutive_previousDotState2(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance)))
           endif
           if (any(numerics_integrator == 4_pInt)) then
             allocate(constitutive_RK4dotState(g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance))) 
           endif
           if (any(numerics_integrator == 5_pInt)) then
             do s = 1_pInt,6_pInt
               allocate(constitutive_RKCK45dotState(s,g,i,e)%p(constitutive_nonlocal_sizeDotState(myInstance))) 
             enddo
           endif
           constitutive_aTolState(g,i,e)%p =        constitutive_nonlocal_aTolState(myInstance)
           constitutive_sizeState(g,i,e) =          constitutive_nonlocal_sizeState(myInstance)
           constitutive_sizeDotState(g,i,e) =       constitutive_nonlocal_sizeDotState(myInstance)
           constitutive_sizePostResults(g,i,e) =    constitutive_nonlocal_sizePostResults(myInstance)
           
         case default
           call IO_error(201_pInt,ext_msg=trim(phase_plasticity(material_phase(g,i,e))))            ! unknown plasticity
          
       end select
     enddo
   enddo
 enddo
 call constitutive_nonlocal_stateInit(constitutive_state0(1,1:iMax,1:eMax))
 do e = 1_pInt,mesh_NcpElems                                                                        ! loop over elements
   myNgrains = homogenization_Ngrains(mesh_element(3,e)) 
   forall(i = 1_pInt:FE_Nips(FE_geomtype(mesh_element(2,e))), g = 1_pInt:myNgrains)
     constitutive_partionedState0(g,i,e)%p = constitutive_state0(g,i,e)%p
     constitutive_state(g,i,e)%p = constitutive_state0(g,i,e)%p                                     ! need to be defined for first call of constitutive_microstructure in crystallite_init
   endforall
 enddo
 
!--------------------------------------------------------------------------------------------------
! write out state size file
 call IO_write_jobBinaryIntFile(777,'sizeStateConst', size(constitutive_sizeState))
 write (777,rec=1) constitutive_sizeState
 close(777)

!--------------------------------------------------------------------------------------------------
! report
 constitutive_maxSizeState       = maxval(constitutive_sizeState)
 constitutive_maxSizeDotState    = maxval(constitutive_sizeDotState)
 constitutive_maxSizePostResults = maxval(constitutive_sizePostResults)
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) then
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_state0:          ', shape(constitutive_state0)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_partionedState0: ', shape(constitutive_partionedState0)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_subState0:       ', shape(constitutive_subState0)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_state:           ', shape(constitutive_state)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_aTolState:       ', shape(constitutive_aTolState)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_dotState:        ', shape(constitutive_dotState)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_deltaState:      ', shape(constitutive_deltaState)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_sizeState:       ', shape(constitutive_sizeState)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_sizeDotState:    ', shape(constitutive_sizeDotState)
   write(6,'(a32,1x,7(i8,1x))') 'constitutive_sizePostResults: ', shape(constitutive_sizePostResults)
   write(6,*)
   write(6,'(a32,1x,7(i8,1x))') 'maxSizeState:       ', constitutive_maxSizeState
   write(6,'(a32,1x,7(i8,1x))') 'maxSizeDotState:    ', constitutive_maxSizeDotState
   write(6,'(a32,1x,7(i8,1x))') 'maxSizePostResults: ', constitutive_maxSizePostResults
 endif
 flush(6)
 
end subroutine constitutive_init


!--------------------------------------------------------------------------------------------------
!> @brief returns the homogenize elasticity matrix
!--------------------------------------------------------------------------------------------------
pure function constitutive_homogenizedC(g,i,e)
 use material, only: &
   phase_plasticity, &
   material_phase
 use constitutive_none
 use constitutive_j2
 use constitutive_phenopowerlaw
 use constitutive_titanmod
 use constitutive_dislotwin
 use constitutive_nonlocal
 
 implicit none
 real(pReal), dimension(6,6) :: constitutive_homogenizedC
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number

 select case (phase_plasticity(material_phase(g,i,e)))
 
   case (constitutive_none_label)
     constitutive_homogenizedC = constitutive_none_homogenizedC(constitutive_state,g,i,e)
     
   case (constitutive_j2_label)
     constitutive_homogenizedC = constitutive_j2_homogenizedC(constitutive_state,g,i,e)
     
   case (constitutive_phenopowerlaw_label)
     constitutive_homogenizedC = constitutive_phenopowerlaw_homogenizedC(constitutive_state,g,i,e)

   case (constitutive_titanmod_label)
     constitutive_homogenizedC = constitutive_titanmod_homogenizedC(constitutive_state,g,i,e)
 
   case (constitutive_dislotwin_label)
     constitutive_homogenizedC = constitutive_dislotwin_homogenizedC(constitutive_state,g,i,e)
     
   case (constitutive_nonlocal_label)
     constitutive_homogenizedC = constitutive_nonlocal_homogenizedC(constitutive_state,g,i,e)
     
 end select

end function constitutive_homogenizedC


!--------------------------------------------------------------------------------------------------
!> @brief returns average length of Burgers vector
!--------------------------------------------------------------------------------------------------
real(pReal) function constitutive_averageBurgers(g,i,e)
 use material, only: &
   phase_plasticity,material_phase
 use constitutive_none, only: &
   constitutive_none_label
 use constitutive_j2, only: &
   constitutive_j2_label
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_label
 use constitutive_titanmod, only: &
   constitutive_titanmod_label
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_label
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_label
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number 
   e                                                                                                !< element number

 select case (phase_plasticity(material_phase(g,i,e)))
 
   case (constitutive_none_label)
     constitutive_averageBurgers = 2.5e-10_pReal                                                    ! constitutive_none_averageBurgers(constitutive_state,g,i,e)
     
   case (constitutive_j2_label)
     constitutive_averageBurgers = 2.5e-10_pReal                                                    ! constitutive_j2_averageBurgers(constitutive_state,g,i,e)
     
   case (constitutive_phenopowerlaw_label)
     constitutive_averageBurgers = 2.5e-10_pReal                                                    ! constitutive_phenopowerlaw_averageBurgers(constitutive_state,g,i,e)

   case (constitutive_titanmod_label)
     constitutive_averageBurgers = 2.5e-10_pReal                                                    ! constitutive_titanmod_averageBurgers(constitutive_state,g,i,e)
     
   case (constitutive_dislotwin_label)
     constitutive_averageBurgers = 2.5e-10_pReal                                                    ! constitutive_dislotwin_averageBurgers(constitutive_state,g,i,e)
     
   case (constitutive_nonlocal_label)
     constitutive_averageBurgers = 2.5e-10_pReal                                                    ! constitutive_nonlocal_averageBurgers(constitutive_state,g,i,e)
     
 end select

end function constitutive_averageBurgers


!--------------------------------------------------------------------------------------------------
!> @brief calls microstructure function of the different constitutive models
!--------------------------------------------------------------------------------------------------
subroutine constitutive_microstructure(Temperature, Fe, Fp, g, i, e)
 use material, only: &
   phase_plasticity, &
   material_phase
 use constitutive_none, only: &
   constitutive_none_label, &
   constitutive_none_microstructure
 use constitutive_j2, only: &
   constitutive_j2_label, &
   constitutive_j2_microstructure
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_label, &
   constitutive_phenopowerlaw_microstructure
 use constitutive_titanmod, only: &
   constitutive_titanmod_label, &
   constitutive_titanmod_microstructure
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_label, &
   constitutive_dislotwin_microstructure
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_label, &
   constitutive_nonlocal_microstructure
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),   intent(in) :: &
   Temperature
 real(pReal),   intent(in), dimension(3,3) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 
 select case (phase_plasticity(material_phase(g,i,e)))
  
   case (constitutive_none_label)
     call constitutive_none_microstructure(Temperature,constitutive_state,g,i,e)
      
   case (constitutive_j2_label)
     call constitutive_j2_microstructure(Temperature,constitutive_state,g,i,e)
      
   case (constitutive_phenopowerlaw_label)
     call constitutive_phenopowerlaw_microstructure(Temperature,constitutive_state,g,i,e)
   
   case (constitutive_titanmod_label)
     call constitutive_titanmod_microstructure(Temperature,constitutive_state,g,i,e)
    
   case (constitutive_dislotwin_label)
     call constitutive_dislotwin_microstructure(Temperature,constitutive_state,g,i,e)
    
   case (constitutive_nonlocal_label)
     call constitutive_nonlocal_microstructure(constitutive_state, Temperature, Fe, Fp, g, i, e)
      
 end select
 
end subroutine constitutive_microstructure


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient                                 *
!--------------------------------------------------------------------------------------------------
subroutine constitutive_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, Temperature, g, i, e)
 use material, only: &
   phase_plasticity, &
   material_phase
 use constitutive_none, only: &
   constitutive_none_label, &
   constitutive_none_LpAndItsTangent
 use constitutive_j2, only: &
   constitutive_j2_label, &
   constitutive_j2_LpAndItsTangent
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_label, &
   constitutive_phenopowerlaw_LpAndItsTangent
 use constitutive_titanmod, only: &
   constitutive_titanmod_label, &
   constitutive_titanmod_LpAndItsTangent
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_label, &
   constitutive_dislotwin_LpAndItsTangent
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_label, &
   constitutive_nonlocal_LpAndItsTangent
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),   intent(in) :: &
   Temperature
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 real(pReal),   intent(out), dimension(3,3) :: &
   Lp                                                                                               !< plastic velocity gradient
 real(pReal),   intent(out), dimension(9,9) :: &
   dLp_dTstar                                                                                       !< derivative of Lp with respect to Tstar (4th-order tensor)
 
 select case (phase_plasticity(material_phase(g,i,e)))
 
   case (constitutive_none_label)
     call constitutive_none_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_j2_label)
     call constitutive_j2_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_phenopowerlaw_label)
     call constitutive_phenopowerlaw_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,g,i,e)
   
   case (constitutive_titanmod_label)
     call constitutive_titanmod_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_dislotwin_label)
     call constitutive_dislotwin_LpAndItsTangent(Lp,dLp_dTstar,Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_nonlocal_label)
     call constitutive_nonlocal_LpAndItsTangent(Lp, dLp_dTstar, Tstar_v, Temperature, constitutive_state(g,i,e), g, i, e)
    
 end select
 
end subroutine constitutive_LpAndItsTangent



!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to 
!> the elastic deformation gradient depending on the selected elastic law
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_TandItsTangent(T, dT_dFe, Fe, g, i, e)
 use material, only: &
   phase_elasticity,material_phase
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe                                                                                               !< elastic deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: &
   dT_dFe                                                                                           !< derivative of 2nd P-K stress with respect to elastic deformation gradient
 
 select case (phase_elasticity(material_phase(g,i,e)))
 
   case (constitutive_hooke_label)
       call constitutive_hooke_TandItsTangent(T, dT_dFe, Fe, g, i, e)
     
 end select
 
end subroutine constitutive_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief returns the 2nd Piola-Kirchhoff stress tensor and its tangent with respect to 
!> the elastic deformation gradient depending for hookes law
!--------------------------------------------------------------------------------------------------
pure subroutine constitutive_hooke_TandItsTangent(T, dT_dFe, Fe, g, i, e)
use math, only : &
  math_mul33x33, &
  math_mul3333xx33, &
  math_Mandel66to3333, &
  math_transpose33, &
  math_I3

 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),   intent(in),  dimension(3,3) :: &
   Fe                                                                                               !< elastic deformation gradient
 real(pReal),   intent(out), dimension(3,3) :: &
   T                                                                                                !< 2nd Piola-Kirchhoff stress tensor
 real(pReal),   intent(out), dimension(3,3,3,3) :: & 
   dT_dFe                                                                                           !< dT/dFe
 
 integer(pInt) :: o, p
 real(pReal), dimension(3,3)     :: FeT
 real(pReal), dimension(3,3,3,3) :: C


 C = math_Mandel66to3333(constitutive_homogenizedC(g,i,e))

 FeT = math_transpose33(Fe)
 T = 0.5_pReal*math_mul3333xx33(C,math_mul33x33(FeT,Fe)-math_I3)

 forall (o=1_pInt:3_pInt, p=1_pInt:3_pInt) dT_dFe(o,p,1:3,1:3) = math_mul33x33(C(o,p,1:3,1:3), FeT) ! dT*_ij/dFe_kl

end subroutine constitutive_hooke_TandItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure 
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDotState(Tstar_v, Fe, Fp, Temperature, subdt, subfrac, g, i, e)
 use prec, only: &
   pLongInt
 use debug, only: &
   debug_cumDotStateCalls, &
   debug_cumDotStateTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips, &
   mesh_maxNipNeighbors
 use material, only: &
   phase_plasticity, &
   material_phase, &
   homogenization_maxNgrains
 use constitutive_none, only: &
   constitutive_none_dotState, &
   constitutive_none_label
 use constitutive_j2, only:  &
   constitutive_j2_dotState, &
   constitutive_j2_label
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_dotState, &
   constitutive_phenopowerlaw_label
 use constitutive_titanmod, only: &
   constitutive_titanmod_dotState, &
   constitutive_titanmod_label
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_dotState, &
   constitutive_dislotwin_label
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_dotState, &
   constitutive_nonlocal_label
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),  intent(in) :: &
   Temperature, &
   subdt                                                                                            !< timestep
 real(pReal),  intent(in), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   subfrac                                                                                          !< subfraction of timestep
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   Fe, &                                                                                            !< elastic deformation gradient
   Fp                                                                                               !< plastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)
 integer(pLongInt) :: &
   tick, tock, & 
   tickrate, &
   maxticks
 
 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 
 select case (phase_plasticity(material_phase(g,i,e)))
 
   case (constitutive_none_label)
     constitutive_dotState(g,i,e)%p = constitutive_none_dotState(Tstar_v,Temperature,constitutive_state,g,i,e)
  
   case (constitutive_j2_label)
     constitutive_dotState(g,i,e)%p = constitutive_j2_dotState(Tstar_v,Temperature,constitutive_state,g,i,e)
  
   case (constitutive_phenopowerlaw_label)
     constitutive_dotState(g,i,e)%p = constitutive_phenopowerlaw_dotState(Tstar_v,Temperature,constitutive_state,g,i,e)
 
   case (constitutive_titanmod_label)
     constitutive_dotState(g,i,e)%p = constitutive_titanmod_dotState(Tstar_v,Temperature,constitutive_state,g,i,e)
   
   case (constitutive_dislotwin_label)
     constitutive_dotState(g,i,e)%p = constitutive_dislotwin_dotState(Tstar_v,Temperature,constitutive_state,g,i,e)
  
   case (constitutive_nonlocal_label)
     constitutive_dotState(g,i,e)%p = constitutive_nonlocal_dotState(Tstar_v, Fe, Fp, Temperature, constitutive_state, &
                                                                       constitutive_state0, subdt, subfrac, g, i, e)
 
 end select
 
 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDotState)
     debug_cumDotStateCalls = debug_cumDotStateCalls + 1_pInt
     debug_cumDotStateTicks = debug_cumDotStateTicks + tock-tick
     !$OMP FLUSH (debug_cumDotStateTicks)
     if (tock < tick) debug_cumDotStateTicks  = debug_cumDotStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDotState)
 endif
 
end subroutine constitutive_collectDotState


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the incremental change of 
!> microstructure based on the current stress and state  
!--------------------------------------------------------------------------------------------------
subroutine constitutive_collectDeltaState(Tstar_v, Temperature, g, i, e)
 use prec, only: &
   pLongInt
 use debug, only: &
   debug_cumDeltaStateCalls, &
   debug_cumDeltaStateTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   phase_plasticity, &
   material_phase, &
   homogenization_maxNgrains
 use constitutive_none, only: &
   constitutive_none_deltaState, &
   constitutive_none_label
 use constitutive_j2, only: &
   constitutive_j2_deltaState, &
   constitutive_j2_label
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_deltaState, &
   constitutive_phenopowerlaw_label
 use constitutive_titanmod, only: &
   constitutive_titanmod_deltaState, &
   constitutive_titanmod_label
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_deltaState, &
   constitutive_dislotwin_label
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_deltaState, &
   constitutive_nonlocal_label
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),   intent(in) :: &
   Temperature
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress     
 integer(pLongInt) :: &
   tick, tock, & 
   tickrate, &
   maxticks

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)

 select case (phase_plasticity(material_phase(g,i,e)))

   case (constitutive_none_label)
     constitutive_deltaState(g,i,e)%p = constitutive_none_deltaState(Tstar_v,Temperature,constitutive_state,g,i,e)
 
   case (constitutive_j2_label)
     constitutive_deltaState(g,i,e)%p = constitutive_j2_deltaState(Tstar_v,Temperature,constitutive_state,g,i,e)
 
   case (constitutive_phenopowerlaw_label)
     constitutive_deltaState(g,i,e)%p = constitutive_phenopowerlaw_deltaState(Tstar_v,Temperature,constitutive_state,g,i,e)

   case (constitutive_titanmod_label)
     constitutive_deltaState(g,i,e)%p = constitutive_titanmod_deltaState(Tstar_v,Temperature,constitutive_state,g,i,e)
  
   case (constitutive_dislotwin_label)
     constitutive_deltaState(g,i,e)%p = constitutive_dislotwin_deltaState(Tstar_v,Temperature,constitutive_state,g,i,e)
 
   case (constitutive_nonlocal_label)
     call constitutive_nonlocal_deltaState(constitutive_deltaState(g,i,e),constitutive_state, Tstar_v,Temperature,g,i,e)
 
 end select

 if (iand(debug_level(debug_constitutive), debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDeltaState)
     debug_cumDeltaStateCalls = debug_cumDeltaStateCalls + 1_pInt
     debug_cumDeltaStateTicks = debug_cumDeltaStateTicks + tock-tick
     !$OMP FLUSH (debug_cumDeltaStateTicks)
     if (tock < tick) debug_cumDeltaStateTicks  = debug_cumDeltaStateTicks + maxticks
   !$OMP END CRITICAL (debugTimingDeltaState)
 endif

end subroutine constitutive_collectDeltaState


!--------------------------------------------------------------------------------------------------
!> @brief contains the constitutive equation for calculating the rate of change of microstructure
!--------------------------------------------------------------------------------------------------
real(pReal) function constitutive_dotTemperature(Tstar_v,Temperature,g,i,e)
 use prec, only: &
   pLongInt
 use debug, only: &
   debug_cumDotTemperatureCalls, &
   debug_cumDotTemperatureTicks, &
   debug_level, &
   debug_constitutive, &
   debug_levelBasic
 use material, only: &
   phase_plasticity, &
   material_phase
 use constitutive_none, only: &
   constitutive_none_dotTemperature, &
   constitutive_none_label
 use constitutive_j2, only: &
   constitutive_j2_dotTemperature, &
   constitutive_j2_label
 use constitutive_phenopowerlaw, only: & 
   constitutive_phenopowerlaw_dotTemperature, &
   constitutive_phenopowerlaw_label
 use constitutive_titanmod, only: &
   constitutive_titanmod_dotTemperature, &
   constitutive_titanmod_label
 use constitutive_dislotwin, only:     constitutive_dislotwin_dotTemperature, &
                                       constitutive_dislotwin_label
 use constitutive_nonlocal, only:      constitutive_nonlocal_dotTemperature, &
                                       constitutive_nonlocal_label

 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal),   intent(in) :: &
   Temperature
 real(pReal),   intent(in),  dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola-Kirchhoff stress
 integer(pLongInt) :: &
   tick, tock, & 
   tickrate, &
   maxticks
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) &
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 
 select case (phase_plasticity(material_phase(g,i,e)))
 
   case (constitutive_none_label)
     constitutive_dotTemperature = constitutive_none_dotTemperature(Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_j2_label)
     constitutive_dotTemperature = constitutive_j2_dotTemperature(Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_phenopowerlaw_label)
     constitutive_dotTemperature = constitutive_phenopowerlaw_dotTemperature(Tstar_v,Temperature,constitutive_state,g,i,e)
 
   case (constitutive_titanmod_label)
     constitutive_dotTemperature = constitutive_titanmod_dotTemperature(Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_dislotwin_label)
     constitutive_dotTemperature = constitutive_dislotwin_dotTemperature(Tstar_v,Temperature,constitutive_state,g,i,e)
    
   case (constitutive_nonlocal_label)
     constitutive_dotTemperature = constitutive_nonlocal_dotTemperature(Tstar_v,Temperature,constitutive_state,g,i,e)
    
 end select
 
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0_pInt) then
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   !$OMP CRITICAL (debugTimingDotTemperature)
     debug_cumDotTemperatureCalls = debug_cumDotTemperatureCalls + 1_pInt
     debug_cumDotTemperatureTicks = debug_cumDotTemperatureTicks + tock-tick
     !$OMP FLUSH (debug_cumDotTemperatureTicks)
     if (tock < tick) debug_cumDotTemperatureTicks  = debug_cumDotTemperatureTicks + maxticks
   !$OMP END CRITICAL (debugTimingDotTemperature)
 endif
 
end function constitutive_dotTemperature


!--------------------------------------------------------------------------------------------------
!> @brief returns array of constitutive results
!--------------------------------------------------------------------------------------------------
function constitutive_postResults(Tstar_v, Fe, Temperature, dt, g, i, e)
 use mesh, only: &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   phase_plasticity, &
   material_phase, &
   homogenization_maxNgrains
 use constitutive_none, only: &
   constitutive_none_postResults, &
   constitutive_none_label
 use constitutive_j2, only: &
   constitutive_j2_postResults, &
   constitutive_j2_label
 use constitutive_phenopowerlaw, only: &
   constitutive_phenopowerlaw_postResults, &
   constitutive_phenopowerlaw_label
 use constitutive_titanmod, only: &
   constitutive_titanmod_postResults, &
   constitutive_titanmod_label
 use constitutive_dislotwin, only: &
   constitutive_dislotwin_postResults, &
   constitutive_dislotwin_label
 use constitutive_nonlocal, only: &
   constitutive_nonlocal_postResults, &
   constitutive_nonlocal_label
 
 implicit none
 integer(pInt), intent(in) :: &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number
 real(pReal), dimension(constitutive_sizePostResults(g,i,e)) :: &
   constitutive_postResults
 real(pReal),  intent(in) :: &
   Temperature, &
   dt                                                                                               !< timestep
 real(pReal),  intent(in), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
   Fe                                                                                               !< elastic deformation gradient
 real(pReal),  intent(in), dimension(6) :: &
   Tstar_v                                                                                          !< 2nd Piola Kirchhoff stress tensor (Mandel)

 constitutive_postResults = 0.0_pReal
 
 select case (phase_plasticity(material_phase(g,i,e)))
 
   case (constitutive_none_label)
     constitutive_postResults = constitutive_none_postResults(Tstar_v,Temperature,dt,constitutive_state,g,i,e)
    
   case (constitutive_j2_label)
     constitutive_postResults = constitutive_j2_postResults(Tstar_v,Temperature,dt,constitutive_state,g,i,e)
    
   case (constitutive_phenopowerlaw_label)
     constitutive_postResults = constitutive_phenopowerlaw_postResults(Tstar_v,Temperature,dt,constitutive_state,g,i,e)
   
   case (constitutive_titanmod_label)
     constitutive_postResults = constitutive_titanmod_postResults(Tstar_v,Temperature,dt,constitutive_state,g,i,e)
    
   case (constitutive_dislotwin_label)
     constitutive_postResults = constitutive_dislotwin_postResults(Tstar_v,Temperature,dt,constitutive_state,g,i,e)
    
   case (constitutive_nonlocal_label)
     constitutive_postResults = constitutive_nonlocal_postResults(Tstar_v, Fe, Temperature, dt, constitutive_state, &
                                                                  constitutive_dotstate(g,i,e), g, i, e)
 end select
  
end function constitutive_postResults


end module constitutive