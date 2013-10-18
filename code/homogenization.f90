! Copyright 2011-13 Max-Planck-Institut für Eisenforschung GmbH
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
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @brief homogenization manager, organizing deformation partitioning and stress homogenization 
!--------------------------------------------------------------------------------------------------
module homogenization
 use prec, only: &
   pInt, &
   pReal, &
   p_vec
 
!--------------------------------------------------------------------------------------------------
! General variables for the homogenization at a  material point
 implicit none
 private
 type(p_vec),   dimension(:,:),         allocatable, public :: &
   homogenization_state0                                                                            !< pointer array to homogenization state at start of FE increment
 real(pReal),   dimension(:,:,:,:),     allocatable, public :: &
   materialpoint_F0, &                                                                              !< def grad of IP at start of FE increment
   materialpoint_F, &                                                                               !< def grad of IP to be reached at end of FE increment
   materialpoint_P                                                                                  !< first P--K stress of IP
 real(pReal),   dimension(:,:,:,:,:,:), allocatable, public ::  &
   materialpoint_dPdF                                                                               !< tangent of first P--K stress at IP
 real(pReal),   dimension(:,:,:),       allocatable, public :: &
   materialpoint_results                                                                            !< results array of material point

 type(p_vec),   dimension(:,:),         allocatable, public, protected :: &
   homogenization_state                                                                             !< pointer array to current homogenization state (end of converged time step)
 integer(pInt), dimension(:,:),         allocatable, public, protected  :: &
   homogenization_sizeState                                                                         !< size of state array per grain
 integer(pInt),                                      public, protected  :: &
   materialpoint_sizeResults, &
   homogenization_maxSizePostResults

 type(p_vec),   dimension(:,:),         allocatable, private :: &
   homogenization_subState0                                                                         !< pointer array to homogenization state at start of homogenization increment
 real(pReal),   dimension(:,:,:,:),     allocatable, private :: &
   materialpoint_subF0, &                                                                           !< def grad of IP at beginning of homogenization increment
   materialpoint_subF                                                                               !< def grad of IP to be reached at end of homog inc
 real(pReal),   dimension(:,:),         allocatable, private :: &
   materialpoint_subFrac, &
   materialpoint_subStep, &
   materialpoint_subdt, &
   materialpoint_heat
 integer(pInt), dimension(:,:),         allocatable, private :: &
   homogenization_sizePostResults                                                                   !< size of postResults array per material point
 integer(pInt),                                      private :: &
   homogenization_maxSizeState
 logical,       dimension(:,:),         allocatable, private :: &
   materialpoint_requested, &
   materialpoint_converged
 logical,       dimension(:,:,:),       allocatable, private :: &
   materialpoint_doneAndHappy

 public ::  &
   homogenization_init, &
   materialpoint_stressAndItsTangent, &
   materialpoint_postResults
 private :: &
   homogenization_partitionDeformation, &
   homogenization_updateState, &
   homogenization_averageStressAndItsTangent, &
   homogenization_averageHeat, &
   homogenization_postResults

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!--------------------------------------------------------------------------------------------------
subroutine homogenization_init()
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use math, only: &
   math_I3
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelBasic, &
   debug_e, &
   debug_g
 use IO, only: &
   IO_error, &
   IO_open_file, &
   IO_open_jobFile_stat, &
   IO_write_jobFile, &
   IO_write_jobIntFile, &
   IO_timeStamp 
 use mesh, only: &
   mesh_maxNips, &
   mesh_NcpElems, &
   mesh_element, & 
   FE_Nips, &
   FE_geomtype
 use constitutive, only: &
   constitutive_maxSizePostResults
 use crystallite, only: &
   crystallite_maxSizePostResults
 use material
 use homogenization_isostrain
 use homogenization_RGC

 implicit none
 integer(pInt), parameter :: fileunit = 200_pInt
 integer(pInt) :: e,i,p,myInstance
 integer(pInt), dimension(:,:), pointer :: thisSize
 character(len=64), dimension(:,:), pointer :: thisOutput
 logical :: knownHomogenization
 
!--------------------------------------------------------------------------------------------------
! parse homogenization from config file
 if (.not. IO_open_jobFile_stat(fileunit,material_localFileExt)) &                                  ! no local material configuration present...
   call IO_open_file(fileunit,material_configFile)                                                  ! ... open material.config file
 call homogenization_isostrain_init(fileunit)
 call homogenization_RGC_init(fileunit)
 close(fileunit)
 
!--------------------------------------------------------------------------------------------------
! write description file for homogenization output
 call IO_write_jobFile(fileunit,'outputHomogenization')
 do p = 1,material_Nhomogenization
   i = homogenization_typeInstance(p)                                                               ! which instance of this homogenization type
   knownHomogenization = .true.                                                                     ! assume valid
   select case(homogenization_type(p))                                                              ! split per homogenization type
     case (homogenization_isostrain_label)
       thisOutput => homogenization_isostrain_output
       thisSize   => homogenization_isostrain_sizePostResult
     case (homogenization_RGC_label)
       thisOutput => homogenization_RGC_output
       thisSize   => homogenization_RGC_sizePostResult
     case default
       knownHomogenization = .false.
   end select   
   write(fileunit,'(/,a,/)')  '['//trim(homogenization_name(p))//']'
   if (knownHomogenization) then
     write(fileunit,'(a)') '(type)'//char(9)//trim(homogenization_type(p))
     write(fileunit,'(a,i4)') '(ngrains)'//char(9),homogenization_Ngrains(p)
     do e = 1,homogenization_Noutput(p)
       write(fileunit,'(a,i4)') trim(thisOutput(e,i))//char(9),thisSize(e,i)
     enddo
   endif  
 enddo
 close(fileunit)
 
!--------------------------------------------------------------------------------------------------
! allocate and initialize global variables
 allocate(homogenization_state0(mesh_maxNips,mesh_NcpElems))
 allocate(homogenization_subState0(mesh_maxNips,mesh_NcpElems))
 allocate(homogenization_state(mesh_maxNips,mesh_NcpElems))
 allocate(homogenization_sizeState(mesh_maxNips,mesh_NcpElems))
          homogenization_sizeState                                    = 0_pInt
 allocate(homogenization_sizePostResults(mesh_maxNips,mesh_NcpElems))
          homogenization_sizePostResults                              = 0_pInt
 allocate(materialpoint_heat(mesh_maxNips,mesh_NcpElems))
          materialpoint_heat                                          = 0.0_pReal
 allocate(materialpoint_dPdF(3,3,3,3,mesh_maxNips,mesh_NcpElems))
          materialpoint_dPdF                                          = 0.0_pReal
 allocate(materialpoint_F0(3,3,mesh_maxNips,mesh_NcpElems))
 allocate(materialpoint_F(3,3,mesh_maxNips,mesh_NcpElems))
          materialpoint_F                                             = 0.0_pReal
 allocate(materialpoint_subF0(3,3,mesh_maxNips,mesh_NcpElems))
          materialpoint_subF0                                         = 0.0_pReal
 allocate(materialpoint_subF(3,3,mesh_maxNips,mesh_NcpElems))
          materialpoint_subF                                          = 0.0_pReal
 allocate(materialpoint_P(3,3,mesh_maxNips,mesh_NcpElems))
          materialpoint_P                                             = 0.0_pReal
 allocate(materialpoint_subFrac(mesh_maxNips,mesh_NcpElems))
          materialpoint_subFrac                                       = 0.0_pReal
 allocate(materialpoint_subStep(mesh_maxNips,mesh_NcpElems))
          materialpoint_subStep                                       = 0.0_pReal
 allocate(materialpoint_subdt(mesh_maxNips,mesh_NcpElems))
          materialpoint_subdt                                         = 0.0_pReal
 allocate(materialpoint_requested(mesh_maxNips,mesh_NcpElems))
          materialpoint_requested                                     = .false.
 allocate(materialpoint_converged(mesh_maxNips,mesh_NcpElems))
          materialpoint_converged                                     = .true.
 allocate(materialpoint_doneAndHappy(2,mesh_maxNips,mesh_NcpElems))
          materialpoint_doneAndHappy                                  = .true.

 materialpoint_F0 = spread(spread(math_I3,3,mesh_maxNips),4,mesh_NcpElems)                          ! initialize to identity
 materialpoint_F  = materialpoint_F0
 
!--------------------------------------------------------------------------------------------------
! allocate and initialize global state and postresutls variables
 elementLooping: do e = 1,mesh_NcpElems
   myInstance = homogenization_typeInstance(mesh_element(3,e))
   IpLooping: do i = 1,FE_Nips(FE_geomtype(mesh_element(2,e)))
     select case(homogenization_type(mesh_element(3,e)))
       case (homogenization_isostrain_label)
         if (homogenization_isostrain_sizeState(myInstance) > 0_pInt) then
           allocate(homogenization_state0(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           allocate(homogenization_subState0(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           allocate(homogenization_state(i,e)%p(homogenization_isostrain_sizeState(myInstance)))
           homogenization_state0(i,e)%p  = 0.0_pReal
           homogenization_sizeState(i,e) = homogenization_isostrain_sizeState(myInstance)
         endif
         homogenization_sizePostResults(i,e) = homogenization_isostrain_sizePostResults(myInstance)
       case (homogenization_RGC_label)
         if (homogenization_RGC_sizeState(myInstance) > 0_pInt) then
           allocate(homogenization_state0(i,e)%p(homogenization_RGC_sizeState(myInstance)))
           allocate(homogenization_subState0(i,e)%p(homogenization_RGC_sizeState(myInstance)))
           allocate(homogenization_state(i,e)%p(homogenization_RGC_sizeState(myInstance)))
           homogenization_state0(i,e)%p  = 0.0_pReal
           homogenization_sizeState(i,e) = homogenization_RGC_sizeState(myInstance)
         endif
         homogenization_sizePostResults(i,e) = homogenization_RGC_sizePostResults(myInstance)
       case default
         call IO_error(500_pInt,ext_msg=homogenization_type(mesh_element(3,e)))                     ! unknown homogenization
     end select
   enddo IpLooping
 enddo elementLooping

!--------------------------------------------------------------------------------------------------
! write state size file out
 call IO_write_jobIntFile(777,'sizeStateHomog',size(homogenization_sizeState))
 write (777,rec=1) homogenization_sizeState
 close(777)
 
 homogenization_maxSizeState       = maxval(homogenization_sizeState)
 homogenization_maxSizePostResults = maxval(homogenization_sizePostResults)  
 materialpoint_sizeResults = 1 &                                                                    ! grain count
                           + 1 + homogenization_maxSizePostResults &                                ! homogSize & homogResult
                           + homogenization_maxNgrains * (1 + crystallite_maxSizePostResults &      ! crystallite size & crystallite results
                                                        + 1 + constitutive_maxSizePostResults)      ! constitutive size & constitutive results
 allocate(materialpoint_results(materialpoint_sizeResults,mesh_maxNips,mesh_NcpElems))
 
 write(6,'(/,a)')   ' <<<+-  homogenization init  -+>>>'
 write(6,'(a)')     ' $Id$'
 write(6,'(a15,a)') ' Current time: ',IO_timeStamp()
#include "compilation_info.f90"
 if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_state0:          ', shape(homogenization_state0)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_subState0:       ', shape(homogenization_subState0)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_state:           ', shape(homogenization_state)
   write(6,'(a32,1x,7(i8,1x))')   'homogenization_sizeState:       ', shape(homogenization_sizeState)
   write(6,'(a32,1x,7(i8,1x),/)') 'homogenization_sizePostResults: ', shape(homogenization_sizePostResults)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_dPdF:             ', shape(materialpoint_dPdF)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F0:               ', shape(materialpoint_F0)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_F:                ', shape(materialpoint_F)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF0:            ', shape(materialpoint_subF0)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subF:             ', shape(materialpoint_subF)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_P:                ', shape(materialpoint_P)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_heat:             ', shape(materialpoint_heat)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subFrac:          ', shape(materialpoint_subFrac)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subStep:          ', shape(materialpoint_subStep)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_subdt:            ', shape(materialpoint_subdt)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_requested:        ', shape(materialpoint_requested)
   write(6,'(a32,1x,7(i8,1x))')   'materialpoint_converged:        ', shape(materialpoint_converged)
   write(6,'(a32,1x,7(i8,1x),/)') 'materialpoint_doneAndHappy:     ', shape(materialpoint_doneAndHappy)
   write(6,'(a32,1x,7(i8,1x),/)') 'materialpoint_results:          ', shape(materialpoint_results)
   write(6,'(a32,1x,7(i8,1x))')   'maxSizeState:       ', homogenization_maxSizeState
   write(6,'(a32,1x,7(i8,1x))')   'maxSizePostResults: ', homogenization_maxSizePostResults
 endif
 flush(6)
 
 if (debug_g < 1 .or. debug_g > homogenization_Ngrains(mesh_element(3,debug_e))) &
   call IO_error(602_pInt,ext_msg='component (grain)')
 
end subroutine homogenization_init


!--------------------------------------------------------------------------------------------------
!> @brief  parallelized calculation of stress and corresponding tangent at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_stressAndItsTangent(updateJaco,dt)
 use numerics, only: &
   subStepMinHomog, &
   subStepSizeHomog, &
   stepIncreaseHomog, &
   nHomog, &
   nMPstate
 use math, only: & 
   math_transpose33
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP, &
   terminallyIll
 use mesh, only: &
   mesh_element, &
   mesh_NcpElems, &
   mesh_maxNips
 use material, only: &
   homogenization_Ngrains
 use constitutive, only: &
   constitutive_state0, &
   constitutive_partionedState0, &
   constitutive_state
 use crystallite, only: &
   crystallite_heat, &
   crystallite_F0, &
   crystallite_Fp0, &
   crystallite_Fp, &
   crystallite_Lp0, &
   crystallite_Lp, &
   crystallite_dPdF, &
   crystallite_dPdF0, &
   crystallite_Tstar0_v, &
   crystallite_Tstar_v, &
   crystallite_partionedF0, &
   crystallite_partionedF, &
   crystallite_partionedFp0, &
   crystallite_partionedLp0, &
   crystallite_partioneddPdF0, &
   crystallite_partionedTstar0_v, &
   crystallite_dt, &
   crystallite_requested, &
   crystallite_converged, &
   crystallite_stressAndItsTangent, &
   crystallite_orientations
 use debug, only: &
   debug_level, &
   debug_homogenization, &
   debug_levelBasic, &
   debug_levelSelective, &
   debug_e, &
   debug_i, &
   debug_MaterialpointLoopDistribution, &
   debug_MaterialpointStateLoopDistribution
 use math, only: &
   math_pDecomposition
                          
 implicit none
 real(pReal), intent(in) :: dt                                                                      !< time increment
 logical,     intent(in) :: updateJaco                                                              !< initiating Jacobian update
 logical                 :: rate_sensitivity
 integer(pInt) :: &
   NiterationHomog, &
   NiterationMPstate, &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e, &                                                                                             !< element number
   myNgrains

!--------------------------------------------------------------------------------------------------
! initialize to starting condition
 if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
     write(6,'(/a,i5,1x,i2)') '<< HOMOG >> Material Point start at el ip ', debug_e, debug_i

     write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F0', &
                                     math_transpose33(materialpoint_F0(1:3,1:3,debug_i,debug_e))
     write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< HOMOG >> F', &
                                     math_transpose33(materialpoint_F(1:3,1:3,debug_i,debug_e))
   !$OMP END CRITICAL (write2out)
 endif

!--------------------------------------------------------------------------------------------------
! initialize restoration points of ...
 do e = FEsolving_execElem(1),FEsolving_execElem(2)
   myNgrains = homogenization_Ngrains(mesh_element(3,e))
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)
     constitutive_partionedState0(g,i,e)%p = constitutive_state0(g,i,e)%p                           ! ...microstructures
     crystallite_partionedFp0(1:3,1:3,g,i,e) = crystallite_Fp0(1:3,1:3,g,i,e)                       ! ...plastic def grads
     crystallite_partionedLp0(1:3,1:3,g,i,e) = crystallite_Lp0(1:3,1:3,g,i,e)                       ! ...plastic velocity grads
     crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,g,i,e) = crystallite_dPdF0(1:3,1:3,1:3,1:3,g,i,e)   ! ...stiffness
     crystallite_partionedF0(1:3,1:3,g,i,e) = crystallite_F0(1:3,1:3,g,i,e)                         ! ...def grads
     crystallite_partionedTstar0_v(1:6,g,i,e) = crystallite_Tstar0_v(1:6,g,i,e)                     ! ...2nd PK stress
   endforall
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e))
     materialpoint_subF0(1:3,1:3,i,e) = materialpoint_F0(1:3,1:3,i,e)                               ! ...def grad
     materialpoint_subFrac(i,e) = 0.0_pReal
     materialpoint_subStep(i,e) = 1.0_pReal/subStepSizeHomog                                        ! <<added to adopt flexibility in cutback size>>
     materialpoint_converged(i,e) = .false.                                                         ! pretend failed step of twice the required size
     materialpoint_requested(i,e) = .true.                                                          ! everybody requires calculation
   endforall
   forall(i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), homogenization_sizeState(i,e) > 0_pInt) &
     homogenization_subState0(i,e)%p = homogenization_state0(i,e)%p                                 ! ...internal homogenization state
 enddo

 NiterationHomog = 0_pInt
 
 cutBackLooping: do while (.not. terminallyIll .and. &
      any(materialpoint_subStep(:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinHomog))

   !$OMP PARALLEL DO PRIVATE(myNgrains)
   elementLooping1: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     IpLooping1: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
              
       converged: if ( materialpoint_converged(i,e) ) then
#ifndef _OPENMP
         if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt &
            .and. ((e == debug_e .and. i == debug_i) & 
                   .or. .not. iand(debug_level(debug_homogenization),debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,1x,f12.8,1x,a,1x,f12.8,1x,a,i8,1x,i2/)') '<< HOMOG >> winding forward from', &
             materialpoint_subFrac(i,e), 'to current materialpoint_subFrac', &
             materialpoint_subFrac(i,e)+materialpoint_subStep(i,e),'in materialpoint_stressAndItsTangent at el ip',e,i
         endif
#endif
         
!--------------------------------------------------------------------------------------------------
! calculate new subStep and new subFrac
         materialpoint_subFrac(i,e) = materialpoint_subFrac(i,e) + materialpoint_subStep(i,e)
         !$OMP FLUSH(materialpoint_subFrac)
         materialpoint_subStep(i,e) = min(1.0_pReal-materialpoint_subFrac(i,e), &
                                          stepIncreaseHomog*materialpoint_subStep(i,e))                   ! introduce flexibility for step increase/acceleration
         !$OMP FLUSH(materialpoint_subStep)
                  
         steppingNeeded: if (materialpoint_subStep(i,e) > subStepMinHomog) then
         
           ! wind forward grain starting point of...
           crystallite_partionedF0(1:3,1:3,1:myNgrains,i,e) = crystallite_partionedF(1:3,1:3,1:myNgrains,i,e) ! ...def grads
           crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e) = crystallite_Fp(1:3,1:3,1:myNgrains,i,e)    ! ...plastic def grads
           crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e) = crystallite_Lp(1:3,1:3,1:myNgrains,i,e)    ! ...plastic velocity grads
           crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,1:myNgrains,i,e) = crystallite_dPdF(1:3,1:3,1:3,1:3,1:myNgrains,i,e)! ...stiffness
           crystallite_partionedTstar0_v(1:6,1:myNgrains,i,e) = crystallite_Tstar_v(1:6,1:myNgrains,i,e)  ! ...2nd PK stress
           forall (g = 1:myNgrains) constitutive_partionedState0(g,i,e)%p = constitutive_state(g,i,e)%p   ! ...microstructures
           if (homogenization_sizeState(i,e) > 0_pInt) &
             homogenization_subState0(i,e)%p = homogenization_state(i,e)%p                                ! ...internal state of homog scheme
           materialpoint_subF0(1:3,1:3,i,e) = materialpoint_subF(1:3,1:3,i,e)                             ! ...def grad
           !$OMP FLUSH(materialpoint_subF0)
         elseif (materialpoint_requested(i,e)) then steppingNeeded                                        ! already at final time (??)
           if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionHomog)
               debug_MaterialpointLoopDistribution(min(nHomog+1,NiterationHomog)) = &
                 debug_MaterialpointLoopDistribution(min(nHomog+1,NiterationHomog)) + 1
             !$OMP END CRITICAL (distributionHomog)
           endif
         endif steppingNeeded

       else converged
         if ( (myNgrains == 1_pInt .and. materialpoint_subStep(i,e) <= 1.0 ) .or. &                         ! single grain already tried internal subStepping in crystallite
              subStepSizeHomog * materialpoint_subStep(i,e) <=  subStepMinHomog ) then                      ! would require too small subStep
                                                                                                            ! cutback makes no sense
           !$OMP FLUSH(terminallyIll)
           if (.not. terminallyIll) then                                                                    ! so first signals terminally ill...
             !$OMP CRITICAL (write2out)
               write(6,*) 'Integration point ', i,' at element ', e, ' terminally ill'
             !$OMP END CRITICAL (write2out)
           endif
           !$OMP CRITICAL (setTerminallyIll)
             terminallyIll = .true.                                                                         ! ...and kills all others
           !$OMP END CRITICAL (setTerminallyIll)
         else                                                                                               ! cutback makes sense
           materialpoint_subStep(i,e) = subStepSizeHomog * materialpoint_subStep(i,e)                       ! crystallite had severe trouble, so do a significant cutback
           !$OMP FLUSH(materialpoint_subStep)
           
#ifndef _OPENMP
           if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt &
              .and. ((e == debug_e .and. i == debug_i) &
                    .or. .not. iand(debug_level(debug_homogenization), debug_levelSelective) /= 0_pInt)) then
             write(6,'(a,1x,f12.8,a,i8,1x,i2/)') &
               '<< HOMOG >> cutback step in materialpoint_stressAndItsTangent with new materialpoint_subStep:',&
               materialpoint_subStep(i,e),' at el ip',e,i 
           endif
#endif
  
!--------------------------------------------------------------------------------------------------
! restore...
           crystallite_Fp(1:3,1:3,1:myNgrains,i,e) = crystallite_partionedFp0(1:3,1:3,1:myNgrains,i,e)      ! ...plastic def grads
           crystallite_Lp(1:3,1:3,1:myNgrains,i,e) = crystallite_partionedLp0(1:3,1:3,1:myNgrains,i,e)      ! ...plastic velocity grads
           crystallite_dPdF(1:3,1:3,1:3,1:3,1:myNgrains,i,e) = crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,1:myNgrains,i,e) ! ...stiffness
           crystallite_Tstar_v(1:6,1:myNgrains,i,e) = crystallite_partionedTstar0_v(1:6,1:myNgrains,i,e)    ! ...2nd PK stress
           forall (g = 1:myNgrains) constitutive_state(g,i,e)%p = constitutive_partionedState0(g,i,e)%p     ! ...microstructures
           if (homogenization_sizeState(i,e) > 0_pInt) &
             homogenization_state(i,e)%p = homogenization_subState0(i,e)%p                                  ! ...internal state of homog scheme
         endif       
       endif converged
     
       if (materialpoint_subStep(i,e) > subStepMinHomog) then
         materialpoint_requested(i,e) = .true.
         materialpoint_subF(1:3,1:3,i,e) = materialpoint_subF0(1:3,1:3,i,e) + &
                                         materialpoint_subStep(i,e) * (materialpoint_F(1:3,1:3,i,e) - materialpoint_F0(1:3,1:3,i,e))
         materialpoint_subdt(i,e) = materialpoint_subStep(i,e) * dt
         materialpoint_doneAndHappy(1:2,i,e) = [.false.,.true.]
       endif
     enddo IpLooping1
   enddo elementLooping1
   !$OMP END PARALLEL DO

   NiterationMPstate = 0_pInt
   
   convergenceLooping: do while (.not. terminallyIll .and. &
             any(            materialpoint_requested(:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                 .and. .not. materialpoint_doneAndHappy(1,:,FEsolving_execELem(1):FEsolving_execElem(2)) &
                ) .and. &
             NiterationMPstate < nMPstate)
     NiterationMPstate = NiterationMPstate + 1

!--------------------------------------------------------------------------------------------------
! deformation partitioning
! based on materialpoint_subF0,.._subF,crystallite_partionedF0, and homogenization_state, 
! results in crystallite_partionedF
     !$OMP PARALLEL DO PRIVATE(myNgrains)
     elementLooping2: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       myNgrains = homogenization_Ngrains(mesh_element(3,e))
       IpLooping2: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if (      materialpoint_requested(i,e) .and. &                                             ! process requested but...
             .not. materialpoint_doneAndHappy(1,i,e)) then                                          ! ...not yet done material points
           call homogenization_partitionDeformation(i,e)                                            ! partition deformation onto constituents
           crystallite_dt(1:myNgrains,i,e) = materialpoint_subdt(i,e)                               ! propagate materialpoint dt to grains
           crystallite_requested(1:myNgrains,i,e) = .true.                                          ! request calculation for constituents
         else
           crystallite_requested(1:myNgrains,i,e) = .false.                                         ! calculation for constituents not required anymore
         endif
       enddo IpLooping2
     enddo elementLooping2
     !$OMP END PARALLEL DO
 
!--------------------------------------------------------------------------------------------------
! crystallite integration
! based on crystallite_partionedF0,.._partionedF
! incrementing by crystallite_dt
     rate_sensitivity = .false.                                                                      ! request rate sensitive contribution to dPdF
     call crystallite_stressAndItsTangent(updateJaco,rate_sensitivity)                               ! request stress and tangent calculation for constituent grains

!--------------------------------------------------------------------------------------------------
! state update
     !$OMP PARALLEL DO
     elementLooping3: do e = FEsolving_execElem(1),FEsolving_execElem(2)
       IpLooping3: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
         if (      materialpoint_requested(i,e) .and. &
             .not. materialpoint_doneAndHappy(1,i,e)) then
           if (.not. all(crystallite_converged(:,i,e))) then
             materialpoint_doneAndHappy(1:2,i,e) = [.true.,.false.]
             materialpoint_converged(i,e) = .false.
           else
             materialpoint_doneAndHappy(1:2,i,e) = homogenization_updateState(i,e)
             materialpoint_converged(i,e) = all(homogenization_updateState(i,e))                     ! converged if done and happy
           endif
           !$OMP FLUSH(materialpoint_converged)
           if (materialpoint_converged(i,e)) then
             if (iand(debug_level(debug_homogenization), debug_levelBasic) /= 0_pInt) then
               !$OMP CRITICAL (distributionMPState)
                 debug_MaterialpointStateLoopdistribution(NiterationMPstate) = &
                   debug_MaterialpointStateLoopdistribution(NiterationMPstate) + 1_pInt
               !$OMP END CRITICAL (distributionMPState)
             endif
           endif
         endif
       enddo IpLooping3
     enddo elementLooping3
     !$OMP END PARALLEL DO

   enddo convergenceLooping

   NiterationHomog = NiterationHomog + 1_pInt

 enddo cutBackLooping

 if (.not. terminallyIll ) then   
   call crystallite_orientations()                                                                  ! calculate crystal orientations
   !$OMP PARALLEL DO
   elementLooping4: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     IpLooping4: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       call homogenization_averageStressAndItsTangent(i,e)
       materialpoint_heat(i,e) = homogenization_averageHeat(i,e)   
     enddo IpLooping4
   enddo elementLooping4
   !$OMP END PARALLEL DO
 else
   !$OMP CRITICAL (write2out)
   write(6,'(/,a,/)') '<< HOMOG >> Material Point terminally ill'
   !$OMP END CRITICAL (write2out)
 endif
 
end subroutine materialpoint_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief parallelized calculation of result array at material points
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_postResults
 use FEsolving, only: &
   FEsolving_execElem, &
   FEsolving_execIP
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_Ngrains, &
   microstructure_crystallite
 use constitutive, only: &
   constitutive_sizePostResults, &
   constitutive_postResults
 use crystallite, only: &
   crystallite_sizePostResults, &
   crystallite_postResults

 implicit none
 integer(pInt) :: &
   thePos, &
   theSize, &
   myNgrains, &
   myCrystallite, &
   g, &                                                                                             !< grain number
   i, &                                                                                             !< integration point number
   e                                                                                                !< element number

 !$OMP PARALLEL DO PRIVATE(myNgrains,myCrystallite,thePos,theSize)
   elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
     myNgrains = homogenization_Ngrains(mesh_element(3,e))
     myCrystallite = microstructure_crystallite(mesh_element(4,e))
     IpLooping: do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       thePos = 0_pInt
       
       theSize = homogenization_sizePostResults(i,e)
       materialpoint_results(thePos+1,i,e) = real(theSize,pReal)                                    ! tell size of homogenization results
       thePos = thePos + 1_pInt

       if (theSize > 0_pInt) then                                                                   ! any homogenization results to mention?
         materialpoint_results(thePos+1:thePos+theSize,i,e) = homogenization_postResults(i,e)       ! tell homogenization results
         thePos = thePos + theSize
       endif
       
       materialpoint_results(thePos+1,i,e) = real(myNgrains,pReal)                                  ! tell number of grains at materialpoint
       thePos = thePos + 1_pInt

       grainLooping :do g = 1,myNgrains
         theSize = (1 + crystallite_sizePostResults(myCrystallite)) + (1 + constitutive_sizePostResults(g,i,e))
         materialpoint_results(thePos+1:thePos+theSize,i,e) = crystallite_postResults(g,i,e)       ! tell crystallite results
         thePos = thePos + theSize
       enddo grainLooping
     enddo IpLooping
   enddo elementLooping
 !$OMP END PARALLEL DO

end subroutine materialpoint_postResults
 
 
!--------------------------------------------------------------------------------------------------
!> @brief  partition material point def grad onto constituents
!--------------------------------------------------------------------------------------------------
subroutine homogenization_partitionDeformation(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains
 use crystallite, only: &
   crystallite_partionedF0, &
   crystallite_partionedF
 use homogenization_isostrain, only: &
   homogenization_isostrain_label, &
   homogenization_isostrain_partitionDeformation
 use homogenization_RGC, only: &
   homogenization_RGC_label, &
   homogenization_RGC_partitionDeformation

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 
 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label) chosenHomogenization
     call homogenization_isostrain_partitionDeformation(&
                          crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                          materialpoint_subF(1:3,1:3,ip,el),&
                          el)
   case (homogenization_RGC_label) chosenHomogenization
     call homogenization_RGC_partitionDeformation(&
                         crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                         materialpoint_subF(1:3,1:3,ip,el),&
                         homogenization_state(ip,el), &
                         ip, &
                         el)
 end select chosenHomogenization

end subroutine homogenization_partitionDeformation


!--------------------------------------------------------------------------------------------------
!> @brief update the internal state of the homogenization scheme and tell whether "done" and 
!> "happy" with result
!--------------------------------------------------------------------------------------------------
function homogenization_updateState(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains
 use crystallite, only: &
   crystallite_P, &
   crystallite_dPdF, &
   crystallite_partionedF,&
   crystallite_partionedF0
 use homogenization_RGC, only: &
   homogenization_RGC_updateState, &
   homogenization_RGC_label
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 logical, dimension(2) :: homogenization_updateState
 
 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))

   case (homogenization_RGC_label) chosenHomogenization
     homogenization_updateState = &
        homogenization_RGC_updateState( homogenization_state(ip,el), &
                                        homogenization_subState0(ip,el), &
                                        crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                        crystallite_partionedF(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                        crystallite_partionedF0(1:3,1:3,1:homogenization_maxNgrains,ip,el),&
                                        materialpoint_subF(1:3,1:3,ip,el),&
                                        materialpoint_subdt(ip,el), &
                                        crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
                                        ip, &
                                        el)
   case default chosenHomogenization
     homogenization_updateState = .true.
 end select chosenHomogenization

end function homogenization_updateState


!--------------------------------------------------------------------------------------------------
!> @brief derive average stress and stiffness from constituent quantities
!--------------------------------------------------------------------------------------------------
subroutine homogenization_averageStressAndItsTangent(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type, &
   homogenization_maxNgrains
 use crystallite, only: &
   crystallite_P,crystallite_dPdF
 use homogenization_isostrain, only: &
   homogenization_isostrain_averageStressAndItsTangent, &
   homogenization_isostrain_label
 use homogenization_RGC, only: &
   homogenization_RGC_averageStressAndItsTangent, &
   homogenization_RGC_label

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 
 chosenHomogenization: select case(homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label) chosenHomogenization
     call homogenization_isostrain_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       el)
   case (homogenization_RGC_label) chosenHomogenization
     call homogenization_RGC_averageStressAndItsTangent(&
       materialpoint_P(1:3,1:3,ip,el), &
       materialpoint_dPdF(1:3,1:3,1:3,1:3,ip,el),&
       crystallite_P(1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       crystallite_dPdF(1:3,1:3,1:3,1:3,1:homogenization_maxNgrains,ip,el), &
       el)
 end select chosenHomogenization

end subroutine homogenization_averageStressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief derive average heat from constituent quantities (does not depend on choosen 
!! homogenization scheme)
!--------------------------------------------------------------------------------------------------
real(pReal) function homogenization_averageHeat(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_Ngrains
 use crystallite, only: &
   crystallite_heat

 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point number
   el                                                                                               !< element number
 integer(pInt) :: &
   Ngrains

!--------------------------------------------------------------------------------------------------
! computing the average heat
 Ngrains = homogenization_Ngrains(mesh_element(3,el))
 homogenization_averageHeat= sum(crystallite_heat(1:Ngrains,ip,el))/real(Ngrains,pReal)

end function homogenization_averageHeat


!--------------------------------------------------------------------------------------------------
!> @brief return array of homogenization results for post file inclusion. call only, 
!> if homogenization_sizePostResults(i,e) > 0 !!
!--------------------------------------------------------------------------------------------------
function homogenization_postResults(ip,el)
 use mesh, only: &
   mesh_element
 use material, only: &
   homogenization_type
 use homogenization_isostrain, only: &
   homogenization_isostrain_postResults, &
   homogenization_isostrain_label
 use homogenization_RGC, only: &
   homogenization_RGC_postResults, &
   homogenization_RGC_label
 
 implicit none
 integer(pInt), intent(in) :: &
   ip, &                                                                                            !< integration point
   el                                                                                               !< element number
 real(pReal), dimension(homogenization_sizePostResults(ip,el)) :: homogenization_postResults

 homogenization_postResults = 0.0_pReal
 chosenHomogenization: select case (homogenization_type(mesh_element(3,el)))
   case (homogenization_isostrain_label) chosenHomogenization
     homogenization_postResults = homogenization_isostrain_postResults(&
                                  ip, &
                                  el, &
                                  materialpoint_P(1:3,1:3,ip,el), &
                                  materialpoint_F(1:3,1:3,ip,el))
   case (homogenization_RGC_label) chosenHomogenization
     homogenization_postResults = homogenization_RGC_postResults(&
                                  homogenization_state(ip,el),&
                                  ip, &
                                  el, &
                                  materialpoint_P(1:3,1:3,ip,el), &
                                  materialpoint_F(1:3,1:3,ip,el))
 end select chosenHomogenization

end function homogenization_postResults

end module homogenization
