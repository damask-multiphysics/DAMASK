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
! $Id$
!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @brief crystallite state integration functions and reporting of results
!--------------------------------------------------------------------------------------------------

module crystallite
 use prec, only: pReal, pInt

 implicit none
 private
 character(len=64), dimension(:,:),        allocatable, private :: &
   crystallite_output                                                                               !< name of each post result output
 integer(pInt),                                         public, protected :: &
   crystallite_maxSizePostResults
 integer(pInt),     dimension(:),          allocatable, public, protected :: &
   crystallite_sizePostResults
 integer(pInt),     dimension(:,:),        allocatable, private :: &
   crystallite_sizePostResult
 integer(pInt),     dimension(:,:,:),      allocatable, private :: &
   crystallite_symmetryID                                                                           !< crystallographic symmetry 1=cubic 2=hexagonal, needed in all orientation calcs     
     
 real(pReal),       dimension (:,:,:),     allocatable, public :: &
   crystallite_dt, &                                                                                !< requested time increment of each grain
   crystallite_Temperature, &                                                                       !< Temp of each grain
   crystallite_partionedTemperature0                                                                !< Temp of each grain at start of homog inc
 real(pReal),       dimension (:,:,:),     allocatable, private :: &
   crystallite_subdt, &                                                                             !< substepped time increment of each grain
   crystallite_subFrac, &                                                                           !< already calculated fraction of increment
   crystallite_subStep, &                                                                           !< size of next integration step
   crystallite_subTemperature0, &                                                                   !< Temp of each grain at start of crystallite inc
   crystallite_dotTemperature                                                                       !< evolution of Temperature of each grain
 real(pReal),       dimension (:,:,:,:),   allocatable, public :: &
   crystallite_Tstar_v, &                                                                           !< current 2nd Piola-Kirchhoff stress vector (end of converged time step)
   crystallite_Tstar0_v, &                                                                          !< 2nd Piola-Kirchhoff stress vector at start of FE inc
   crystallite_partionedTstar0_v                                                                    !< 2nd Piola-Kirchhoff stress vector at start of homog inc
 real(pReal),       dimension (:,:,:,:),   allocatable, private :: &
   crystallite_subTstar0_v, &                                                                       !< 2nd Piola-Kirchhoff stress vector at start of crystallite inc
   crystallite_orientation, &                                                                       !< orientation as quaternion
   crystallite_orientation0, &                                                                      !< initial orientation as quaternion
   crystallite_rotation                                                                             !< grain rotation away from initial orientation as axis-angle (in degrees) in crystal reference frame 
 real(pReal),       dimension (:,:,:,:,:), allocatable, public :: &
   crystallite_Fp, &                                                                                !< current plastic def grad (end of converged time step)
   crystallite_Fp0, &                                                                               !< plastic def grad at start of FE inc
   crystallite_partionedFp0,&                                                                       !< plastic def grad at start of homog inc
   crystallite_F0, &                                                                                !< def grad at start of FE inc
   crystallite_partionedF,  &                                                                       !< def grad to be reached at end of homog inc
   crystallite_partionedF0, &                                                                       !< def grad at start of homog inc
   crystallite_Lp, &                                                                                !< current plastic velocitiy grad (end of converged time step)
   crystallite_Lp0, &                                                                               !< plastic velocitiy grad at start of FE inc
   crystallite_partionedLp0,&                                                                       !< plastic velocity grad at start of homog inc
   crystallite_P                                                                                    !< 1st Piola-Kirchhoff stress per grain
 real(pReal), dimension (:,:,:,:,:),       allocatable, private :: &
   crystallite_Fe, &                                                                                !< current "elastic" def grad (end of converged time step)
   crystallite_subFe0,&                                                                             !< "elastic" def grad at start of crystallite inc
   crystallite_invFp, &                                                                             !< inverse of current plastic def grad (end of converged time step)
   crystallite_subFp0,&                                                                             !< plastic def grad at start of crystallite inc
   crystallite_subF,  &                                                                             !< def grad to be reached at end of crystallite inc
   crystallite_subF0, &                                                                             !< def grad at start of crystallite inc
   crystallite_subLp0,&                                                                             !< plastic velocity grad at start of crystallite inc
   crystallite_disorientation                                                                       !< disorientation between two neighboring ips (only calculated for single grain IPs)
 real(pReal), dimension (:,:,:,:,:,:,:),   allocatable, public :: &
   crystallite_dPdF, &                                                                              !< current individual dPdF per grain (end of converged time step)
   crystallite_dPdF0, &                                                                             !< individual dPdF per grain at start of FE inc
   crystallite_partioneddPdF0                                                                       !< individual dPdF per grain at start of homog inc
 real(pReal), dimension (:,:,:,:,:,:,:),   allocatable, private :: &
   crystallite_fallbackdPdF                                                                         !< dPdF fallback for non-converged grains (elastic prediction)
 logical,     dimension (:,:,:),           allocatable, public :: &
   crystallite_requested                                                                            !< flag to request crystallite calculation
 logical,     dimension (:,:,:),           allocatable, public, protected :: &
   crystallite_converged                                                                            !< convergence flag
 logical,     dimension (:,:,:),           allocatable, private :: &
   crystallite_localPlasticity, &                                                                   !< indicates this grain to have purely local constitutive law
   crystallite_todo                                                                                 !< flag to indicate need for further computation
 logical,     dimension (:,:),             allocatable, private :: &
   crystallite_clearToWindForward, &
   crystallite_clearToCutback, &
   crystallite_syncSubFrac, &
   crystallite_syncSubFracCompleted, &
   crystallite_neighborEnforcedCutback

 public :: &
   crystallite_init, &
   crystallite_stressAndItsTangent, &
   crystallite_orientations, &
   crystallite_postResults
 private :: &
   crystallite_integrateStateFPI, &   
   crystallite_integrateStateEuler, &
   crystallite_integrateStateAdaptiveEuler, &
   crystallite_integrateStateRK4, &
   crystallite_integrateStateRKCK45, &
   crystallite_integrateStress, &
   crystallite_stateJump
  
contains


!--------------------------------------------------------------------------------------------------
!> @brief allocates and initialize per grain variables
!--------------------------------------------------------------------------------------------------
subroutine crystallite_init(Temperature)
 use, intrinsic :: iso_fortran_env                                                                  ! to get compiler_version and compiler_options (at least for gfortran 4.6 at the moment)
 use debug, only:      debug_info, &
                       debug_reset, &
                       debug_level, &
                       debug_crystallite, &
                       debug_levelBasic
 use math, only:       math_I3, &
                       math_EulerToR, &
                       math_inv33, &
                       math_transpose33, &
                       math_mul33xx33, &
                       math_mul33x33
 use FEsolving, only:  FEsolving_execElem, &
                       FEsolving_execIP
 use mesh, only:       mesh_element, &
                       mesh_NcpElems, &
                       mesh_maxNips, &
                       mesh_maxNipNeighbors
 use IO
 use material
 use lattice, only:    lattice_symmetryType
 
 use constitutive, only: constitutive_microstructure
 use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_label, &
                                       constitutive_phenopowerlaw_structureName
 use constitutive_titanmod, only:      constitutive_titanmod_label, &
                                       constitutive_titanmod_structureName
 use constitutive_dislotwin, only:     constitutive_dislotwin_label, &
                                       constitutive_dislotwin_structureName
 use constitutive_nonlocal, only:      constitutive_nonlocal_label, &
                                       constitutive_nonlocal_structureName
  
 implicit none
 real(pReal), intent(in)  :: Temperature
 integer(pInt), parameter :: myFile = 200_pInt, &
                             maxNchunks = 2_pInt
 
 !*** local variables ***!
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt)               g, &                                                                    ! grain number
                             i, &                                                                    ! integration point number
                             e, &                                                                    ! element number
                             gMax, &                                                                 ! maximum number of grains
                             iMax, &                                                                 ! maximum number of integration points
                             eMax, &                                                                 ! maximum number of elements
                             nMax, &                                                                 ! maximum number of ip neighbors
                             myNgrains, &                                                            ! number of grains in current IP
                             section, &
                             j, &
                             p, &
                             output, &
                             mySize, &
                             myPhase, &
                             myMat
 character(len=64)           tag
 character(len=1024)         line
 
 write(6,'(/,a)') ' <<<+-  crystallite init  -+>>>'
 write(6,'(a)')   ' $Id$'
 write(6,'(a16,a)')   ' Current time : ',IO_timeStamp()
#include "compilation_info.f90"
 
 
 gMax = homogenization_maxNgrains
 iMax = mesh_maxNips
 eMax = mesh_NcpElems
 nMax = mesh_maxNipNeighbors
 
 allocate(crystallite_Temperature(gMax,iMax,eMax));                     crystallite_Temperature = Temperature
 allocate(crystallite_partionedTemperature0(gMax,iMax,eMax)); crystallite_partionedTemperature0 = 0.0_pReal
 allocate(crystallite_subTemperature0(gMax,iMax,eMax));             crystallite_subTemperature0 = 0.0_pReal
 allocate(crystallite_dotTemperature(gMax,iMax,eMax));               crystallite_dotTemperature = 0.0_pReal
 allocate(crystallite_Tstar0_v(6,gMax,iMax,eMax));                         crystallite_Tstar0_v = 0.0_pReal
 allocate(crystallite_partionedTstar0_v(6,gMax,iMax,eMax));       crystallite_partionedTstar0_v = 0.0_pReal
 allocate(crystallite_subTstar0_v(6,gMax,iMax,eMax));                   crystallite_subTstar0_v = 0.0_pReal
 allocate(crystallite_Tstar_v(6,gMax,iMax,eMax));                           crystallite_Tstar_v = 0.0_pReal
 allocate(crystallite_P(3,3,gMax,iMax,eMax));                                     crystallite_P = 0.0_pReal
 allocate(crystallite_F0(3,3,gMax,iMax,eMax));                                   crystallite_F0 = 0.0_pReal
 allocate(crystallite_partionedF0(3,3,gMax,iMax,eMax));                 crystallite_partionedF0 = 0.0_pReal
 allocate(crystallite_partionedF(3,3,gMax,iMax,eMax));                   crystallite_partionedF = 0.0_pReal
 allocate(crystallite_subF0(3,3,gMax,iMax,eMax));                             crystallite_subF0 = 0.0_pReal
 allocate(crystallite_subF(3,3,gMax,iMax,eMax));                               crystallite_subF = 0.0_pReal
 allocate(crystallite_Fp0(3,3,gMax,iMax,eMax));                                 crystallite_Fp0 = 0.0_pReal
 allocate(crystallite_partionedFp0(3,3,gMax,iMax,eMax));               crystallite_partionedFp0 = 0.0_pReal
 allocate(crystallite_subFp0(3,3,gMax,iMax,eMax));                           crystallite_subFp0 = 0.0_pReal
 allocate(crystallite_Fp(3,3,gMax,iMax,eMax));                                   crystallite_Fp = 0.0_pReal
 allocate(crystallite_invFp(3,3,gMax,iMax,eMax));                             crystallite_invFp = 0.0_pReal
 allocate(crystallite_Fe(3,3,gMax,iMax,eMax));                                   crystallite_Fe = 0.0_pReal
 allocate(crystallite_subFe0(3,3,gMax,iMax,eMax));                           crystallite_subFe0 = 0.0_pReal
 allocate(crystallite_Lp0(3,3,gMax,iMax,eMax));                                 crystallite_Lp0 = 0.0_pReal
 allocate(crystallite_partionedLp0(3,3,gMax,iMax,eMax));               crystallite_partionedLp0 = 0.0_pReal
 allocate(crystallite_subLp0(3,3,gMax,iMax,eMax));                           crystallite_subLp0 = 0.0_pReal
 allocate(crystallite_Lp(3,3,gMax,iMax,eMax));                                   crystallite_Lp = 0.0_pReal
 allocate(crystallite_dPdF(3,3,3,3,gMax,iMax,eMax));                           crystallite_dPdF = 0.0_pReal
 allocate(crystallite_dPdF0(3,3,3,3,gMax,iMax,eMax));                         crystallite_dPdF0 = 0.0_pReal
 allocate(crystallite_partioneddPdF0(3,3,3,3,gMax,iMax,eMax));       crystallite_partioneddPdF0 = 0.0_pReal
 allocate(crystallite_fallbackdPdF(3,3,3,3,gMax,iMax,eMax));           crystallite_fallbackdPdF = 0.0_pReal
 allocate(crystallite_dt(gMax,iMax,eMax));                                       crystallite_dt = 0.0_pReal
 allocate(crystallite_subdt(gMax,iMax,eMax));                                 crystallite_subdt = 0.0_pReal
 allocate(crystallite_subFrac(gMax,iMax,eMax));                             crystallite_subFrac = 0.0_pReal
 allocate(crystallite_subStep(gMax,iMax,eMax));                             crystallite_subStep = 0.0_pReal
 allocate(crystallite_orientation(4,gMax,iMax,eMax));                   crystallite_orientation = 0.0_pReal
 allocate(crystallite_orientation0(4,gMax,iMax,eMax));                 crystallite_orientation0 = 0.0_pReal
 allocate(crystallite_rotation(4,gMax,iMax,eMax));                         crystallite_rotation = 0.0_pReal
 allocate(crystallite_disorientation(4,nMax,gMax,iMax,eMax));        crystallite_disorientation = 0.0_pReal
 allocate(crystallite_symmetryID(gMax,iMax,eMax));                       crystallite_symmetryID = 0_pInt
 allocate(crystallite_localPlasticity(gMax,iMax,eMax));             crystallite_localPlasticity = .true.
 allocate(crystallite_requested(gMax,iMax,eMax));                         crystallite_requested = .false.
 allocate(crystallite_todo(gMax,iMax,eMax));                                   crystallite_todo = .false.
 allocate(crystallite_converged(gMax,iMax,eMax));                         crystallite_converged = .true.
 allocate(crystallite_clearToWindForward(iMax,eMax));            crystallite_clearToWindForward = .true.
 allocate(crystallite_syncSubFrac(iMax,eMax));                          crystallite_syncSubFrac = .false.
 allocate(crystallite_syncSubFracCompleted(iMax,eMax));        crystallite_syncSubFracCompleted = .false.
 allocate(crystallite_clearToCutback(iMax,eMax));                    crystallite_clearToCutback = .true.
 allocate(crystallite_neighborEnforcedCutback(iMax,eMax));  crystallite_neighborEnforcedCutback = .false.
 allocate(crystallite_output(maxval(crystallite_Noutput), &
                             material_Ncrystallite)) ;                       crystallite_output = ''
 allocate(crystallite_sizePostResults(material_Ncrystallite)) ;     crystallite_sizePostResults = 0_pInt
 allocate(crystallite_sizePostResult(maxval(crystallite_Noutput), &
                                     material_Ncrystallite)) ;       crystallite_sizePostResult = 0_pInt
 
 
 if (.not. IO_open_jobFile_stat(myFile,material_localFileExt)) then                                 ! no local material configuration present...
   call IO_open_file(myFile,material_configFile)                                                    ! ...open material.config file
 endif
 line = ''
 section = 0_pInt
 
 do while (IO_lc(IO_getTag(line,'<','>')) /= material_partCrystallite)                              ! wind forward to <crystallite>
   read(myFile,'(a1024)',END=100) line
 enddo
 
 do                                                                                                 ! read thru sections of phase part
   read(myFile,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                                                                      ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                                                          ! stop at next part
   if (IO_getTag(line,'[',']') /= '') then                                                          ! next section
     section = section + 1_pInt
     output = 0_pInt                                                                                ! reset output counter
   endif
   if (section > 0_pInt) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1_pInt))                                             ! extract key
     select case(tag)
       case ('(output)')
         output = output + 1_pInt
         crystallite_output(output,section) = IO_lc(IO_stringValue(line,positions,2_pInt))
     end select
   endif
 enddo
 
100 close(myFile)
 
 do i = 1_pInt,material_Ncrystallite                                                                ! sanity checks
 enddo
 
 do i = 1_pInt,material_Ncrystallite
   do j = 1_pInt,crystallite_Noutput(i)
     select case(crystallite_output(j,i))
       case('phase','texture','volume','grainrotationx','grainrotationy','grainrotationz')
         mySize = 1_pInt
       case('orientation','grainrotation')                                                          ! orientation as quaternion, or deviation from initial grain orientation in axis-angle form (angle in degrees)
         mySize = 4_pInt
       case('eulerangles')                                                                          ! Bunge (3-1-3) Euler angles
         mySize = 3_pInt
       case('defgrad','f','fe','fp','lp','e','ee','p','firstpiola','1stpiola','s','tstar','secondpiola','2ndpiola')
         mySize = 9_pInt
       case('elasmatrix')
         mySize = 36_pInt     
       case default
         mySize = 0_pInt     
     end select
   
     if (mySize > 0_pInt) then                                                                      ! any meaningful output found
       crystallite_sizePostResult(j,i) = mySize
       crystallite_sizePostResults(i) = crystallite_sizePostResults(i) + mySize
     endif
   enddo
 enddo
 
 crystallite_maxSizePostResults = 0_pInt
 do j = 1_pInt,material_Nmicrostructure
   if (microstructure_active(j)) &
     crystallite_maxSizePostResults = max(crystallite_maxSizePostResults,&
                                          crystallite_sizePostResults(microstructure_crystallite(j)))
 enddo

! write description file for crystallite output
 
 call IO_write_jobFile(myFile,'outputCrystallite')
  
 do p = 1_pInt,material_Ncrystallite
   write(myFile,*)
   write(myFile,'(a)') '['//trim(crystallite_name(p))//']'
   write(myFile,*)
   do e = 1_pInt,crystallite_Noutput(p)
     write(myFile,'(a,i4)') trim(crystallite_output(e,p))//char(9),crystallite_sizePostResult(e,p)
   enddo
 enddo
 
 close(myFile)
 

do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                      ! iterate over all cp elements
  myNgrains = homogenization_Ngrains(mesh_element(3,e))                                                 ! look up homogenization-->grainCount
  forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1_pInt:myNgrains)
    crystallite_Fp0(1:3,1:3,g,i,e) = math_EulerToR(material_EulerAngles(1:3,g,i,e))                     ! plastic def gradient reflects init orientation
    crystallite_F0(1:3,1:3,g,i,e)  = math_I3
    crystallite_localPlasticity(g,i,e) = phase_localPlasticity(material_phase(g,i,e))
    crystallite_Fe(1:3,1:3,g,i,e)  = math_transpose33(crystallite_Fp0(1:3,1:3,g,i,e))
    crystallite_Fp(1:3,1:3,g,i,e)  = crystallite_Fp0(1:3,1:3,g,i,e)
    crystallite_requested(g,i,e) = .true.
  endforall
enddo
crystallite_partionedTemperature0 = Temperature   ! isothermal assumption
crystallite_partionedFp0 = crystallite_Fp0
crystallite_partionedF0 = crystallite_F0
crystallite_partionedF = crystallite_F0


! Initialize crystallite_symmetryID(g,i,e)

do e = FEsolving_execElem(1),FEsolving_execElem(2)
  myNgrains = homogenization_Ngrains(mesh_element(3,e))
  do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
    do g = 1_pInt,myNgrains
      myPhase = material_phase(g,i,e)
      myMat   = phase_plasticityInstance(myPhase)
      select case (phase_plasticity(myPhase))
        case (constitutive_phenopowerlaw_label)
          crystallite_symmetryID(g,i,e) = lattice_symmetryType(constitutive_phenopowerlaw_structureName(myMat))
        case (constitutive_titanmod_label)
          crystallite_symmetryID(g,i,e) = lattice_symmetryType(constitutive_titanmod_structureName(myMat))
        case (constitutive_dislotwin_label)
          crystallite_symmetryID(g,i,e) = lattice_symmetryType(constitutive_dislotwin_structureName(myMat))
        case (constitutive_nonlocal_label)
          crystallite_symmetryID(g,i,e) = lattice_symmetryType(constitutive_nonlocal_structureName(myMat))
        case default
          crystallite_symmetryID(g,i,e) = 0_pInt ! does this happen for j2 material?
      end select
    enddo
  enddo
enddo

 
call crystallite_orientations()
crystallite_orientation0 = crystallite_orientation             ! Store initial orientations for calculation of grain rotations

!$OMP PARALLEL DO PRIVATE(myNgrains)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do g = 1_pInt,myNgrains
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                         crystallite_Fp(1:3,1:3,g,i,e), g, i, e)  ! update dependent state variables to be consistent with basic states    
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

call crystallite_stressAndItsTangent(.true.,.false.)                   ! request elastic answers
crystallite_fallbackdPdF = crystallite_dPdF                    ! use initial elastic stiffness as fallback
 
!    *** Output  ***
if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Temperature:           ', shape(crystallite_Temperature)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_dotTemperature:        ', shape(crystallite_dotTemperature)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fe:                    ', shape(crystallite_Fe)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fp:                    ', shape(crystallite_Fp)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Lp:                    ', shape(crystallite_Lp)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_F0:                    ', shape(crystallite_F0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Fp0:                   ', shape(crystallite_Fp0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Lp0:                   ', shape(crystallite_Lp0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedF:            ', shape(crystallite_partionedF)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedTemp0:        ', shape(crystallite_partionedTemperature0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedF0:           ', shape(crystallite_partionedF0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedFp0:          ', shape(crystallite_partionedFp0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedLp0:          ', shape(crystallite_partionedLp0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subF:                  ', shape(crystallite_subF)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subTemperature0:       ', shape(crystallite_subTemperature0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_symmetryID:            ', shape(crystallite_symmetryID)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subF0:                 ', shape(crystallite_subF0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFe0:                ', shape(crystallite_subFe0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFp0:                ', shape(crystallite_subFp0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subLp0:                ', shape(crystallite_subLp0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_P:                     ', shape(crystallite_P)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Tstar_v:               ', shape(crystallite_Tstar_v)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_Tstar0_v:              ', shape(crystallite_Tstar0_v)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partionedTstar0_v:     ', shape(crystallite_partionedTstar0_v)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subTstar0_v:           ', shape(crystallite_subTstar0_v)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_dPdF:                  ', shape(crystallite_dPdF)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_dPdF0:                 ', shape(crystallite_dPdF0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_partioneddPdF0:        ', shape(crystallite_partioneddPdF0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_fallbackdPdF:          ', shape(crystallite_fallbackdPdF)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_orientation:           ', shape(crystallite_orientation)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_orientation0:          ', shape(crystallite_orientation0)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_rotation:              ', shape(crystallite_rotation)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_disorientation:        ', shape(crystallite_disorientation)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_dt:                    ', shape(crystallite_dt)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subdt:                 ', shape(crystallite_subdt)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subFrac:               ', shape(crystallite_subFrac)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_subStep:               ', shape(crystallite_subStep)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_localPlasticity:       ', shape(crystallite_localPlasticity)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_requested:             ', shape(crystallite_requested)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_todo:                  ', shape(crystallite_todo)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_converged:             ', shape(crystallite_converged)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_sizePostResults:       ', shape(crystallite_sizePostResults)
    write(6,'(a35,1x,7(i8,1x))') 'crystallite_sizePostResult:        ', shape(crystallite_sizePostResult)
    write(6,*)
    write(6,*) 'Number of nonlocal grains: ',count(.not. crystallite_localPlasticity)
    flush(6)
endif

 call debug_info
 call debug_reset

end subroutine crystallite_init


!--------------------------------------------------------------------------------------------------
!> @brief calculate stress (P) and tangent (dPdF) for crystallites
!--------------------------------------------------------------------------------------------------
subroutine crystallite_stressAndItsTangent(updateJaco,rate_sensitivity)
use numerics, only:      subStepMinCryst, &
                         subStepSizeCryst, &
                         stepIncreaseCryst, &
                         pert_Fg, &
                         pert_method, &
                         nCryst, &
                         numerics_integrator, &
                         numerics_integrationMode, &
                         numerics_timeSyncing, &
                         relevantStrain, &
                         analyticJaco
use debug, only:         debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_CrystalliteLoopDistribution
use IO, only:            IO_warning
use math, only:          math_inv33, &
                         math_identity2nd, &
                         math_transpose33, &
                         math_mul33x33, &
                         math_mul66x6, &
                         math_Mandel6to33, &
                         math_Mandel33to6, &
                         math_I3, &
                         math_mul3333xx3333
use FEsolving, only:     FEsolving_execElem, & 
                         FEsolving_execIP
use mesh, only:          mesh_element, &
                         mesh_NcpElems, &
                         mesh_maxNips, &
                         mesh_ipNeighborhood, &
                         FE_NipNeighbors, &
                         FE_geomtype
use material, only:      homogenization_Ngrains, &
                         homogenization_maxNgrains
use constitutive, only:  constitutive_sizeState, &
                         constitutive_sizeDotState, &
                         constitutive_state, &
                         constitutive_state_backup, &
                         constitutive_subState0, &
                         constitutive_partionedState0, &
                         constitutive_homogenizedC, &
                         constitutive_dotState, &
                         constitutive_dotState_backup, &
                         constitutive_TandItsTangent


implicit none
logical, intent(in) ::            updateJaco, rate_sensitivity                                      ! flag indicating wehther we want to update the Jacobian (stiffness) or not
real(pReal)                       myPert, &                                                         ! perturbation with correct sign
                                  formerSubStep, &
                                  subFracIntermediate
real(pReal), dimension(3,3) ::    invFp, &                                                          ! inverse of the plastic deformation gradient
                                  Fe_guess, &                                                       ! guess for elastic deformation gradient
                                  Tstar                                                             ! 2nd Piola-Kirchhoff stress tensor
real(pReal), dimension(3,3,3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                  dPdF_perturbation1, &
                                  dPdF_perturbation2
real(pReal), dimension(3,3,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                   F_backup, &
                                   Fp_backup, &
                                   InvFp_backup, &
                                   Fe_backup, &
                                   Lp_backup, &
                                   P_backup
real(pReal), dimension(6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                   Tstar_v_backup
real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                   Temperature_backup
integer(pInt)                      NiterationCrystallite, &                                         ! number of iterations in crystallite loop
                                   e, &                                                             ! element index
                                   i, &                                                             ! integration point index
                                   g, &                                                             ! grain index
                                   k, &
                                   l, &
                                   n, &
                                   neighboring_e, &
                                   neighboring_i, &
                                   o, &
                                   p, &
                                   perturbation , &                                                 ! loop counter for forward,backward perturbation mode
                                   myNgrains
logical, dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                      convergenceFlag_backup
! local variables used for calculating analytic Jacobian
real(pReal), dimension(3,3)::       Fpinv_rate, &
                                    FDot_inv, &
                                    junk
real(pReal), dimension(3,3,3,3) ::   dSdFe, &
                                     dFedF, &
                                     dFedFdot, &
                                     dSdF, &
                                     dSdFdot, &
                                     dFp_invdFdot, &
                                     junk2
real(pReal) ::   counter

! --+>> INITIALIZE TO STARTING CONDITION <<+--

if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt&
                        .and. debug_e > 0 .and. debug_e <= mesh_NcpElems &
                        .and. debug_i > 0 .and. debug_i <= mesh_maxNips &
                        .and. debug_g > 0 .and. debug_g <= homogenization_maxNgrains) then
  !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> crystallite start at el ip g ', debug_e, debug_i, debug_g
    write(6,'(a,/,12x,f14.9)') '<< CRYST >> Temp0', crystallite_partionedTemperature0(debug_g,debug_i,debug_e)
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> F0 ', &
                                          math_transpose33(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fp0', &
                                          math_transpose33(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
    write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Lp0', &
                                          math_transpose33(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
  !$OMP END CRITICAL (write2out)
endif

crystallite_subStep = 0.0_pReal

do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                      ! iterate over elements to be processed
  myNgrains = homogenization_Ngrains(mesh_element(3,e))
  forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1_pInt:myNgrains, crystallite_requested(g,i,e))
    crystallite_subTemperature0(g,i,e) = crystallite_partionedTemperature0(g,i,e)                       ! ...temperature
    constitutive_subState0(g,i,e)%p = constitutive_partionedState0(g,i,e)%p                             ! ...microstructure
    crystallite_subFp0(1:3,1:3,g,i,e) = crystallite_partionedFp0(1:3,1:3,g,i,e)                         ! ...plastic def grad
    crystallite_subLp0(1:3,1:3,g,i,e) = crystallite_partionedLp0(1:3,1:3,g,i,e)                         ! ...plastic velocity grad
    crystallite_dPdF0(1:3,1:3,1:3,1:3,g,i,e) = crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,g,i,e)        ! ...stiffness
    crystallite_subF0(1:3,1:3,g,i,e) = crystallite_partionedF0(1:3,1:3,g,i,e)                           ! ...def grad
    crystallite_subTstar0_v(1:6,g,i,e) = crystallite_partionedTstar0_v(1:6,g,i,e)                       !...2nd PK stress
    crystallite_subFe0(1:3,1:3,g,i,e) = math_mul33x33(crystallite_subF0(1:3,1:3,g,i,e), &
                                                      math_inv33(crystallite_subFp0(1:3,1:3,g,i,e)))    ! only needed later on for stiffness calculation
    crystallite_subFrac(g,i,e) = 0.0_pReal
    crystallite_subStep(g,i,e) = 1.0_pReal/subStepSizeCryst
    crystallite_todo(g,i,e) = .true.
    crystallite_converged(g,i,e) = .false.                                                              ! pretend failed step of twice the required size
  endforall
enddo


! --+>> CRYSTALLITE CUTBACK LOOP <<+--

NiterationCrystallite = 0_pInt
numerics_integrationMode = 1_pInt
do while (any(crystallite_todo(:,:,FEsolving_execELem(1):FEsolving_execElem(2))))                       ! cutback loop for crystallites

  if (any(.not. crystallite_localPlasticity) .and. numerics_timeSyncing) then
    
    ! Time synchronization can only be used for nonlocal calculations, and only there it makes sense.
    ! The idea is that in nonlocal calculations often the vast amjority of the ips
    ! converges in one iteration whereas a small fraction of ips has to do a lot of cutbacks. 
    ! Hence, we try to minimize the computational effort by just doing a lot of cutbacks
    ! in the vicinity of the "bad" ips and leave the easily converged volume more or less as it is.
    ! However, some synchronization of the time step has to be done at the border between "bad" ips 
    ! and the ones that immediately converged. 

    if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
      !$OMP CRITICAL (write2out)
      write(6,'(a,i6)') '<< CRYST >> crystallite iteration ',NiterationCrystallite
      !$OMP END CRITICAL (write2out)
    endif

    if (any(crystallite_syncSubFrac)) then 
      
      ! Just did a time synchronization.
      ! If all synchrnizers converged, then do nothing else than winding them forward.
      ! If any of the cynchronizers did not converge, something went completely wrong 
      ! and its not clear how to fix this, so all nonlocals become terminally ill.
      
      if (any(crystallite_syncSubFrac .and. .not. crystallite_converged(1,:,:))) then
        if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
          do e = FEsolving_execElem(1),FEsolving_execElem(2)
            myNgrains = homogenization_Ngrains(mesh_element(3,e))
            do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
              if (crystallite_syncSubFrac(i,e) .and. .not. crystallite_converged(1,i,e)) then
                !$OMP CRITICAL (write2out)
                write(6,'(a,i8,1x,i2)') '<< CRYST >> time synchronization: failed at el,ip ',e,i
                !$OMP END CRITICAL (write2out)
              endif
            enddo
          enddo
        endif
        crystallite_syncSubFrac = .false.
        where(.not. crystallite_localPlasticity)
          crystallite_substep = 0.0_pReal
          crystallite_todo = .false.
        endwhere
      else
        crystallite_clearToWindForward = crystallite_localPlasticity(1,:,:) .or. crystallite_syncSubFrac
        crystallite_clearToCutback = crystallite_localPlasticity(1,:,:)
        if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
          !$OMP CRITICAL (write2out)
          write(6,'(a,i6)') '<< CRYST >> time synchronization: wind forward'
          !$OMP END CRITICAL (write2out)
        endif
      endif

    elseif (any(crystallite_syncSubFracCompleted)) then 
      
      ! Just completed a time synchronization.
      ! Make sure that the ips that synchronized their time step start non-converged 

      where(crystallite_syncSubFracCompleted) &
        crystallite_converged(1,:,:) = .false.
      crystallite_syncSubFracCompleted = .false.
      crystallite_clearToWindForward = crystallite_localPlasticity(1,:,:)
      crystallite_clearToCutback = crystallite_localPlasticity(1,:,:) .or. .not. crystallite_converged(1,:,:)
      if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
        !$OMP CRITICAL (write2out)
        write(6,'(a,i6)') '<< CRYST >> time synchronization: done, proceed with cutback'
        !$OMP END CRITICAL (write2out)
      endif

    else
      
      ! Normal calculation.
      ! If all converged and are at the end of the time increment, then just do a final wind forward.
      ! If all converged, but not all reached the end of the time increment, then we only wind
      ! those forward that are still on their way, all others have to wait.
      ! If some did not converge and all are still at the start of the time increment, 
      ! then all non-convergers force their converged neighbors to also do a cutback.
      ! In case that some ips have already wound forward to an intermediate time (subfrac), 
      ! then all those ips that converged in the first iteration, but now have a non-converged neighbor
      ! have to synchronize their time step to the same intermediate time. If such a synchronization
      ! takes place, all other ips have to wait and only the synchronizers do a cutback. In the next 
      ! iteration those will do a wind forward while all others still wait. 

      crystallite_clearToWindForward = crystallite_localPlasticity(1,:,:)
      crystallite_clearToCutback = crystallite_localPlasticity(1,:,:)
      if (all(crystallite_localPlasticity .or. crystallite_converged)) then
        if (all(crystallite_localPlasticity .or. crystallite_subStep + crystallite_subFrac >= 1.0_pReal)) then
          crystallite_clearToWindForward = .true.   ! final wind forward
          if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
            !$OMP CRITICAL (write2out)
            write(6,'(a,i6)') '<< CRYST >> final wind forward'
            !$OMP END CRITICAL (write2out)
          endif
        else
          crystallite_clearToWindForward = crystallite_localPlasticity(1,:,:) .or. crystallite_subStep(1,:,:) < 1.0_pReal
          if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
            !$OMP CRITICAL (write2out)
            write(6,'(a,i6)') '<< CRYST >> wind forward'
            !$OMP END CRITICAL (write2out)
          endif
        endif
      else
        subFracIntermediate = maxval(crystallite_subFrac, mask=.not.crystallite_localPlasticity)
        if (subFracIntermediate == 0.0_pReal) then
          crystallite_neighborEnforcedCutback = .false.  ! look for ips that require a cutback because of a nonconverged neighbor
          !$OMP PARALLEL DO PRIVATE(neighboring_e,neighboring_i)
            do e = FEsolving_execElem(1),FEsolving_execElem(2) 
              do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
                if (.not. crystallite_localPlasticity(1,i,e) .and. crystallite_converged(1,i,e)) then
                  do n = 1_pInt,FE_NipNeighbors(FE_geomtype(mesh_element(2,e)))
                    neighboring_e = mesh_ipNeighborhood(1,n,i,e)
                    neighboring_i = mesh_ipNeighborhood(2,n,i,e)
                    if (neighboring_e > 0_pInt .and. neighboring_i > 0_pInt) then
                      if (.not. crystallite_localPlasticity(1,neighboring_i,neighboring_e) &
                          .and. .not. crystallite_converged(1,neighboring_i,neighboring_e)) then 
                        crystallite_neighborEnforcedCutback(i,e) = .true.
#ifndef _OPENMP
                        if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
                          write(6,'(a12,i5,1x,i2,a,i5,1x,i2)') '<< CRYST >> ', neighboring_e,neighboring_i, &
                                                                 ' enforced cutback at ',e,i
#endif
                      endif
                    endif
                  enddo
                endif
              enddo
            enddo
          !$OMP END PARALLEL DO
          where(crystallite_neighborEnforcedCutback) &
            crystallite_converged(1,:,:) = .false.
        else
          crystallite_syncSubFrac = .false.  ! look for ips that have to do a time synchronization because of a nonconverged neighbor
          !$OMP PARALLEL DO PRIVATE(neighboring_e,neighboring_i)
            do e = FEsolving_execElem(1),FEsolving_execElem(2) 
              do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
                if (.not. crystallite_localPlasticity(1,i,e) .and. crystallite_subFrac(1,i,e) == 0.0_pReal) then
                  do n = 1_pInt,FE_NipNeighbors(FE_geomtype(mesh_element(2,e)))
                    neighboring_e = mesh_ipNeighborhood(1,n,i,e)
                    neighboring_i = mesh_ipNeighborhood(2,n,i,e)
                    if (neighboring_e > 0_pInt .and. neighboring_i > 0_pInt) then
                      if (.not. crystallite_localPlasticity(1,neighboring_i,neighboring_e) &
                          .and. .not. crystallite_converged(1,neighboring_i,neighboring_e)) then 
                        crystallite_syncSubFrac(i,e) = .true.
#ifndef _OPENMP
                        if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) &
                          write(6,'(a12,i5,1x,i2,a,i5,1x,i2)') '<< CRYST >> ',neighboring_e,neighboring_i, &
                                                               ' enforced time synchronization at ',e,i
#endif
                      endif
                    endif
                  enddo
                endif
              enddo
            enddo
          !$OMP END PARALLEL DO
          where(crystallite_syncSubFrac) &
            crystallite_converged(1,:,:) = .false.
        endif
        where(.not. crystallite_localPlasticity .and. crystallite_subStep < 1.0_pReal) &
          crystallite_converged = .false.
        if (any(crystallite_syncSubFrac)) then   ! have to do syncing now, so all wait except for the synchronizers which do a cutback
          crystallite_clearToWindForward = crystallite_localPlasticity(1,:,:)
          crystallite_clearToCutback = crystallite_localPlasticity(1,:,:) .or. crystallite_syncSubFrac
          if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
            !$OMP CRITICAL (write2out)
            write(6,'(a,i6)') '<< CRYST >> time synchronization: cutback'
            !$OMP END CRITICAL (write2out)
          endif
        else
          where(.not. crystallite_converged(1,:,:)) &
            crystallite_clearToCutback = .true.
          if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
            !$OMP CRITICAL (write2out)
            write(6,'(a,i6)') '<< CRYST >> cutback'
            !$OMP END CRITICAL (write2out)
          endif
        endif
      endif
    endif

    ! Make sure that all cutbackers start with the same substep

    where(.not. crystallite_localPlasticity .and. .not. crystallite_converged) &
      crystallite_subStep = minval(crystallite_subStep, mask=.not. crystallite_localPlasticity &
                                                             .and. .not. crystallite_converged)

    ! Those that do neither wind forward nor cutback are not to do

    where(.not. crystallite_clearToWindForward .and. .not. crystallite_clearToCutback) &
      crystallite_todo(1,:,:) = .false.
      
  endif

  !$OMP PARALLEL DO PRIVATE(myNgrains,formerSubStep)
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                  ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          
          ! --- wind forward ---
          
          if (crystallite_converged(g,i,e) .and. crystallite_clearToWindForward(i,e)) then
            formerSubStep = crystallite_subStep(g,i,e)
            crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
            !$OMP FLUSH(crystallite_subFrac)
            crystallite_subStep(g,i,e) = min(1.0_pReal - crystallite_subFrac(g,i,e), &
                                             stepIncreaseCryst * crystallite_subStep(g,i,e))
            !$OMP FLUSH(crystallite_subStep)
            if (crystallite_subStep(g,i,e) > 0.0_pReal) then
              crystallite_subTemperature0(g,i,e) = crystallite_Temperature(g,i,e)                       ! wind forward...
              crystallite_subF0(1:3,1:3,g,i,e) = crystallite_subF(1:3,1:3,g,i,e)                        ! ...def grad
              !$OMP FLUSH(crystallite_subF0)
              crystallite_subFp0(1:3,1:3,g,i,e) = crystallite_Fp(1:3,1:3,g,i,e)                         ! ...plastic def grad
              crystallite_subFe0(1:3,1:3,g,i,e) = math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), crystallite_invFp(1:3,1:3,g,i,e)) ! only needed later on for stiffness calculation
              crystallite_subLp0(1:3,1:3,g,i,e) = crystallite_Lp(1:3,1:3,g,i,e)                         ! ...plastic velocity gradient
              constitutive_subState0(g,i,e)%p = constitutive_state(g,i,e)%p                             ! ...microstructure
              crystallite_subTstar0_v(1:6,g,i,e) = crystallite_Tstar_v(1:6,g,i,e)                       ! ...2nd PK stress
              if (crystallite_syncSubFrac(i,e)) then                                                    ! if we just did a synchronization of states, then we wind forward without any further time integration
                crystallite_syncSubFracCompleted(i,e) = .true.
                crystallite_syncSubFrac(i,e) = .false.
                crystallite_todo(g,i,e) = .false.
              else
                crystallite_todo(g,i,e) = .true.
              endif
              !$OMP FLUSH(crystallite_todo)
#ifndef _OPENMP
              if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt &
                  .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                         .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
                write(6,'(a,f12.8,a,f12.8,a,i8,1x,i2,1x,i3)') '<< CRYST >> winding forward from ', &
                  crystallite_subFrac(g,i,e)-formerSubStep,' to current crystallite_subfrac ', &
                  crystallite_subFrac(g,i,e),' in crystallite_stressAndItsTangent at el ip g ',e,i,g
                write(6,*)
              endif
#endif
            elseif (formerSubStep > 0.0_pReal) then                                                     ! this crystallite just converged for the entire timestep
              crystallite_todo(g,i,e) = .false.                                                         ! so done here
              !$OMP FLUSH(crystallite_todo)
              if (iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt) then
                !$OMP CRITICAL (distributionCrystallite)
                  debug_CrystalliteLoopDistribution(min(nCryst+1_pInt,NiterationCrystallite)) = &
                    debug_CrystalliteLoopDistribution(min(nCryst+1_pInt,NiterationCrystallite)) + 1_pInt
                !$OMP END CRITICAL (distributionCrystallite)
              endif
            endif
          
          ! --- cutback ---
          
          elseif (.not. crystallite_converged(g,i,e) .and. crystallite_clearToCutback(i,e)) then
            if (crystallite_syncSubFrac(i,e)) then                                                      ! synchronize time
              crystallite_subStep(g,i,e) = subFracIntermediate
            else
              crystallite_subStep(g,i,e) = subStepSizeCryst * crystallite_subStep(g,i,e)                ! cut step in half and restore...
            endif
            !$OMP FLUSH(crystallite_subStep)
            crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e)                         ! ...temperature
            crystallite_Fp(1:3,1:3,g,i,e) = crystallite_subFp0(1:3,1:3,g,i,e)                           ! ...plastic def grad
            !$OMP FLUSH(crystallite_Fp)
            crystallite_invFp(1:3,1:3,g,i,e) = math_inv33(crystallite_Fp(1:3,1:3,g,i,e))
            !$OMP FLUSH(crystallite_invFp)
            crystallite_Lp(1:3,1:3,g,i,e)    = crystallite_subLp0(1:3,1:3,g,i,e)                      ! ...plastic velocity grad
            constitutive_state(g,i,e)%p      = constitutive_subState0(g,i,e)%p                        ! ...microstructure
            crystallite_Tstar_v(1:6,g,i,e)   = crystallite_subTstar0_v(1:6,g,i,e)                     ! ...2nd PK stress
                                                                                                      ! cant restore dotState here, since not yet calculated in first cutback after initialization
            crystallite_todo(g,i,e) = crystallite_subStep(g,i,e) > subStepMinCryst                    ! still on track or already done (beyond repair)
            !$OMP FLUSH(crystallite_todo)
#ifndef _OPENMP
            if(iand(debug_level(debug_crystallite),debug_levelBasic) /= 0_pInt &
                .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
              if (crystallite_todo(g,i,e)) then
                write(6,'(a,f12.8,a,i8,1x,i2,1x,i3)') '<< CRYST >> cutback step in crystallite_stressAndItsTangent &
                                                       &with new crystallite_subStep: ',&
                                                      crystallite_subStep(g,i,e),' at el ip g ',e,i,g
              else
                write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> reached minimum step size &
                                              &in crystallite_stressAndItsTangent at el ip g ',e,i,g
              endif
              write(6,*)
            endif
#endif
          endif

          ! --- prepare for integration ---
            
          if (crystallite_todo(g,i,e) .and. (crystallite_clearToWindForward(i,e) .or. crystallite_clearToCutback(i,e))) then
            crystallite_subF(1:3,1:3,g,i,e) = crystallite_subF0(1:3,1:3,g,i,e) &
                                            + crystallite_subStep(g,i,e) &
                                              * (crystallite_partionedF(1:3,1:3,g,i,e) - crystallite_partionedF0(1:3,1:3,g,i,e))
            !$OMP FLUSH(crystallite_subF)
            crystallite_Fe(1:3,1:3,g,i,e) = math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), crystallite_invFp(1:3,1:3,g,i,e))
            crystallite_subdt(g,i,e) = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
            crystallite_converged(g,i,e) = .false.                                                    ! start out non-converged
          endif
          
        enddo    ! grains
      enddo      ! IPs
    enddo        ! elements
  !$OMP END PARALLEL DO

  if(numerics_timeSyncing) then
    if (any(.not. crystallite_localPlasticity .and. .not. crystallite_todo .and. .not. crystallite_converged &
            .and. crystallite_subStep <= subStepMinCryst)) then                                          ! no way of rescuing a nonlocal ip that violated the lower time step limit, ...
      if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
        do e = FEsolving_execElem(1),FEsolving_execElem(2)
          myNgrains = homogenization_Ngrains(mesh_element(3,e))
          do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
            do g = 1,myNgrains
              if (.not. crystallite_localPlasticity(g,i,e) .and. .not. crystallite_todo(g,i,e) &
                  .and. .not. crystallite_converged(g,i,e) .and. crystallite_subStep(g,i,e) <= subStepMinCryst) then
                !$OMP CRITICAL (write2out)
                write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> nonlocal violated minimum subStep at el,ip,g ',e,i,g
                !$OMP END CRITICAL (write2out)
              endif
            enddo
          enddo
        enddo
      endif
      where(.not. crystallite_localPlasticity)
        crystallite_todo = .false.                                                                       ! ... so let all nonlocal ips die peacefully
        crystallite_subStep = 0.0_pReal
      endwhere
    endif
  endif

  if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt) then
    !$OMP CRITICAL (write2out)
    write(6,*)
    write(6,'(a,e12.5)') '<< CRYST >> min(subStep) ',minval(crystallite_subStep)
    write(6,'(a,e12.5)') '<< CRYST >> max(subStep) ',maxval(crystallite_subStep)
    write(6,'(a,e12.5)') '<< CRYST >> min(subFrac) ',minval(crystallite_subFrac)
    write(6,'(a,e12.5)') '<< CRYST >> max(subFrac) ',maxval(crystallite_subFrac)
    write(6,*)
    !$OMP END CRITICAL (write2out)
  endif

  ! --- integrate --- requires fully defined state array (basic + dependent state)

  if (any(crystallite_todo)) then
    select case(numerics_integrator(numerics_integrationMode))
      case(1_pInt)
        call crystallite_integrateStateFPI()
      case(2_pInt)
        call crystallite_integrateStateEuler()
      case(3_pInt)
        call crystallite_integrateStateAdaptiveEuler()
      case(4_pInt)
        call crystallite_integrateStateRK4()
      case(5_pInt)
        call crystallite_integrateStateRKCK45()
    end select
  endif
  
  where(.not. crystallite_converged .and. crystallite_subStep > subStepMinCryst) &                    ! do not try non-converged & fully cutbacked any further
    crystallite_todo = .true.

  NiterationCrystallite = NiterationCrystallite + 1_pInt
      
enddo                                                                                                 ! cutback loop


! --+>> CHECK FOR NON-CONVERGED CRYSTALLITES <<+--

do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
  myNgrains = homogenization_Ngrains(mesh_element(3,e))
  do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
    do g = 1,myNgrains
      if (.not. crystallite_converged(g,i,e)) then                                                    ! respond fully elastically (might be not required due to becoming terminally ill anyway)
        if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
          !$OMP CRITICAL (write2out)
          write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> no convergence: respond fully elastic at el ip g ',e,i,g
          write(6,*) 
          !$OMP END CRITICAL (write2out)
        endif
        invFp = math_inv33(crystallite_partionedFp0(1:3,1:3,g,i,e))
        Fe_guess = math_mul33x33(crystallite_partionedF(1:3,1:3,g,i,e), invFp)
        call constitutive_TandItsTangent(Tstar, junk2, Fe_guess,g,i,e)
        crystallite_P(1:3,1:3,g,i,e) = math_mul33x33(Fe_guess,math_mul33x33(Tstar,transpose(invFp)))
      endif
      if(iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt &
          .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                  .or. .not. iand(debug_level(debug_crystallite),debug_levelSelective) /= 0_pInt)) then
        !$OMP CRITICAL (write2out)
        write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> central solution of cryst_StressAndTangent at el ip g ',e,i,g
        write(6,*) 
        write(6,'(a,/,3(12x,3(f12.4,1x)/))') '<< CRYST >> P / MPa', math_transpose33(crystallite_P(1:3,1:3,g,i,e))/1.0e6_pReal
        write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Fp', math_transpose33(crystallite_Fp(1:3,1:3,g,i,e))
        write(6,'(a,/,3(12x,3(f14.9,1x)/))') '<< CRYST >> Lp', math_transpose33(crystallite_Lp(1:3,1:3,g,i,e))
        write(6,*) 
        !$OMP END CRITICAL (write2out)
      endif
    enddo
  enddo
enddo


! --+>> STIFFNESS CALCULATION <<+--

if(updateJaco) then                                                                                     ! Jacobian required
  
  if (.not. analyticJaco) then                                                                          ! Calculate Jacobian using perturbations
  
  numerics_integrationMode = 2_pInt
  
  ! --- BACKUP ---
  
  do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)
      constitutive_state_backup(g,i,e)%p(1:constitutive_sizeState(g,i,e)) = &
        constitutive_state(g,i,e)%p(1:constitutive_sizeState(g,i,e))    ! remember unperturbed, converged state, ...
      constitutive_dotState_backup(g,i,e)%p(1:constitutive_sizeDotState(g,i,e)) = &
        constitutive_dotState(g,i,e)%p(1:constitutive_sizeDotState(g,i,e)) ! ... dotStates, ...
    endforall
  enddo
  Temperature_backup     = crystallite_Temperature                                                      ! ... Temperature, ...
  F_backup               = crystallite_subF                                                             ! ... and kinematics
  Fp_backup              = crystallite_Fp
  InvFp_backup           = crystallite_invFp
  Fe_backup              = crystallite_Fe
  Lp_backup              = crystallite_Lp
  Tstar_v_backup         = crystallite_Tstar_v
  P_backup               = crystallite_P
  convergenceFlag_backup = crystallite_converged

  
  ! --- CALCULATE STATE AND STRESS FOR PERTURBATION ---
  
  dPdF_perturbation1 = crystallite_dPdF0                                                            ! initialize stiffness with known good values from last increment
  dPdF_perturbation2 = crystallite_dPdF0                                                            ! initialize stiffness with known good values from last increment
  do perturbation = 1,2                                                                             ! forward and backward perturbation
    if (iand(pert_method,perturbation) > 0_pInt) then                                               ! mask for desired direction
      myPert = -pert_Fg * (-1.0_pReal)**perturbation                                                ! set perturbation step
      do k = 1,3; do l = 1,3                                                                        ! ...alter individual components
        if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
          !$OMP CRITICAL (write2out)
            write(6,'(a,2(1x,i1),1x,a)') '<< CRYST >> [[[[[[ Stiffness perturbation',k,l,']]]]]]'
            write(6,*)
          !$OMP END CRITICAL (write2out)
        endif
        

        ! --- INITIALIZE UNPERTURBED STATE ---
            
        select case(numerics_integrator(numerics_integrationMode))
          case(1_pInt)     ! Fix-point method: restore to last converged state at end of subinc, since this is probably closest to perturbed state
            do e = FEsolving_execElem(1),FEsolving_execElem(2)
              myNgrains = homogenization_Ngrains(mesh_element(3,e))
              forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)
                constitutive_state(g,i,e)%p(1:constitutive_sizeState(g,i,e)) = &
                  constitutive_state_backup(g,i,e)%p(1:constitutive_sizeState(g,i,e))
                constitutive_dotState(g,i,e)%p(1:constitutive_sizeDotState(g,i,e)) = &
                  constitutive_dotState_backup(g,i,e)%p(1:constitutive_sizeDotState(g,i,e))
              endforall
            enddo
            crystallite_Temperature = Temperature_backup
            crystallite_Fp          = Fp_backup 
            crystallite_invFp       = InvFp_backup
            crystallite_Fe          = Fe_backup
            crystallite_Lp          = Lp_backup
            crystallite_Tstar_v     = Tstar_v_backup
          case(2_pInt,3_pInt)                                                                       ! explicit Euler methods: nothing to restore (except for F), since we are only doing a stress integration step
          case(4_pInt,5_pInt)                                                                       ! explicit Runge-Kutta methods: restore to start of subinc, since we are doing a full integration of state and stress
            do e = FEsolving_execElem(1),FEsolving_execElem(2)
              myNgrains = homogenization_Ngrains(mesh_element(3,e))
              forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)
                constitutive_state(g,i,e)%p(1:constitutive_sizeState(g,i,e)) = &
                  constitutive_subState0(g,i,e)%p(1:constitutive_sizeState(g,i,e))
                constitutive_dotState(g,i,e)%p(1:constitutive_sizeDotState(g,i,e)) = & 
                  constitutive_dotState_backup(g,i,e)%p(1:constitutive_sizeDotState(g,i,e))
              endforall
            enddo
            crystallite_Temperature = crystallite_subTemperature0
            crystallite_Fp          = crystallite_subFp0
            crystallite_Fe          = crystallite_subFe0
            crystallite_Lp          = crystallite_subLp0
            crystallite_Tstar_v     = crystallite_subTstar0_v
        end select


        ! --- PERTURB EITHER FORWARD OR BACKWARD ---

        crystallite_subF = F_backup
        crystallite_subF(k,l,:,:,:) = crystallite_subF(k,l,:,:,:) + myPert
        
        crystallite_converged = convergenceFlag_backup
        crystallite_todo = crystallite_requested .and. crystallite_converged
        where (crystallite_todo) crystallite_converged = .false.                                    ! start out non-converged

        select case(numerics_integrator(numerics_integrationMode))
          case(1_pInt)
            call crystallite_integrateStateFPI()
          case(2_pInt)
            call crystallite_integrateStateEuler()
          case(3_pInt)
            call crystallite_integrateStateAdaptiveEuler()
          case(4_pInt)
            call crystallite_integrateStateRK4()
          case(5_pInt)
            call crystallite_integrateStateRKCK45()
        end select
        
        do e = FEsolving_execElem(1),FEsolving_execElem(2)
          myNgrains = homogenization_Ngrains(mesh_element(3,e))
          select case(perturbation)
            case(1_pInt)
              forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                      crystallite_requested(g,i,e) .and. crystallite_converged(g,i,e)) &            ! converged state warrants stiffness update
                dPdF_perturbation1(1:3,1:3,k,l,g,i,e) = (crystallite_P(1:3,1:3,g,i,e) - P_backup(1:3,1:3,g,i,e)) / myPert ! tangent dP_ij/dFg_kl
            case(2_pInt)
              forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                      crystallite_requested(g,i,e) .and. crystallite_converged(g,i,e)) &            ! converged state warrants stiffness update
                dPdF_perturbation2(1:3,1:3,k,l,g,i,e) = (crystallite_P(1:3,1:3,g,i,e) - P_backup(1:3,1:3,g,i,e)) / myPert ! tangent dP_ij/dFg_kl
            end select
        enddo
        
      enddo; enddo                                                                                  ! k,l component perturbation loop

    endif
  enddo                                                                                                 ! perturbation direction


  ! --- STIFFNESS ACCORDING TO PERTURBATION METHOD AND CONVERGENCE ---

  elementLooping: do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    select case(pert_method)
      case(1_pInt)
        forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                crystallite_requested(g,i,e) .and. convergenceFlag_backup(g,i,e)) &                 ! perturbation mode 1: central solution converged
          crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = dPdF_perturbation1(1:3,1:3,1:3,1:3,g,i,e)
      case(2_pInt)
        forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                crystallite_requested(g,i,e) .and. convergenceFlag_backup(g,i,e)) &                 ! perturbation mode 2: central solution converged
          crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = dPdF_perturbation2(1:3,1:3,1:3,1:3,g,i,e)
      case(3_pInt)
        forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
                crystallite_requested(g,i,e) .and. convergenceFlag_backup(g,i,e)) &                 ! perturbation mode 3: central solution converged
          crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = 0.5_pReal* (  dPdF_perturbation1(1:3,1:3,1:3,1:3,g,i,e) &
                                                                + dPdF_perturbation2(1:3,1:3,1:3,1:3,g,i,e))
    end select
    forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains, &
            crystallite_requested(g,i,e) .and. .not. convergenceFlag_backup(g,i,e)) &               ! for any pertubation mode: if central solution did not converge...
      crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = crystallite_fallbackdPdF(1:3,1:3,1:3,1:3,g,i,e)     ! ...use (elastic) fallback
  enddo elementLooping
  

  ! --- RESTORE ---
  
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    forall (i = FEsolving_execIP(1,e):FEsolving_execIP(2,e), g = 1:myNgrains)
      constitutive_state(g,i,e)%p(1:constitutive_sizeState(g,i,e)) = &
        constitutive_state_backup(g,i,e)%p(1:constitutive_sizeState(g,i,e))
      constitutive_dotState(g,i,e)%p(1:constitutive_sizeDotState(g,i,e)) = &
        constitutive_dotState_backup(g,i,e)%p(1:constitutive_sizeDotState(g,i,e))
    endforall
  enddo
  crystallite_Temperature = Temperature_backup
  crystallite_subF        = F_backup
  crystallite_Fp          = Fp_backup 
  crystallite_invFp       = InvFp_backup
  crystallite_Fe          = Fe_backup
  crystallite_Lp          = Lp_backup
  crystallite_Tstar_v     = Tstar_v_backup
  crystallite_P           = P_backup
  crystallite_converged   = convergenceFlag_backup
  
  else                                                                                              ! Calculate Jacobian using analytical expression

  ! --- CALCULATE  ANALYTIC dPdF ---
  
  !$OMP PARALLEL DO PRIVATE(dFedF,dSdF,dSdFe,myNgrains)
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                              ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                            ! iterate over IPs of this element to be processed
        do g = 1_pInt,myNgrains
          dFedF = 0.0_pReal
          do p=1_pInt,3_pInt; do o=1_pInt,3_pInt
            dFedF(p,o,o,1:3) = crystallite_invFp(1:3,p,g,i,e)                                       ! dFe^T_ij/dF_kl = delta_jk * (Fp current^-1)_li
          enddo; enddo
          call constitutive_TandItsTangent(junk,dSdFe,crystallite_subFe0(1:3,1:3,g,i,e),g,i,e)      ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative
          dSdF = math_mul3333xx3333(dSdFe,dFedF)                                                    ! dS/dF = dS/dFe * dFe/dF
          do p=1_pInt,3_pInt; do o=1_pInt,3_pInt
            crystallite_dPdF(1:3,1:3,o,p,g,i,e) = math_mul33x33(math_mul33x33(dFedF(1:3,1:3,o,p),&
                   math_Mandel6to33(crystallite_Tstar_v)),math_transpose33(&
                   crystallite_invFp(1:3,1:3,g,i,e))) &                                             ! dP/dF = dFe/dF * S * Fp^-T...
                   + math_mul33x33(crystallite_subFe0(1:3,1:3,g,i,e),&
                   math_mul33x33(dSdF(1:3,1:3,o,p),math_transpose33(crystallite_invFp(1:3,1:3,g,i,e)))) !         + Fe * dS/dF * Fp^-T         
          enddo; enddo
    enddo; enddo; enddo 
  !$OMP END PARALLEL DO
  endif
  
  if (rate_sensitivity) then
  !$OMP PARALLEL DO PRIVATE(dFedFdot,dSdFdot,dSdFe,Fpinv_rate,FDot_inv,counter,dFp_invdFdot,myNgrains)
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                  ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                ! iterate over IPs of this element to be processed
        do g = 1_pInt,myNgrains
          Fpinv_rate = math_mul33x33(crystallite_invFp(1:3,1:3,g,i,e),crystallite_Lp(1:3,1:3,g,i,e))    ! dFp^-1 = dFp^-1/dt *dt... dFp may overshoot dF by small ammount as 
          FDot_inv = crystallite_subF(1:3,1:3,g,i,e) - crystallite_F0(1:3,1:3,g,i,e)
          counter = 0.0_pReal
          do p=1_pInt,3_pInt; do o=1_pInt,3_pInt
            if (abs(FDot_inv(o,p)) < relevantStrain) then
              FDot_inv(o,p) = 0.0_pReal
            else  
              counter = counter + 1.0_pReal
              FDot_inv(o,p) = crystallite_dt(g,i,e)/FDot_inv(o,p)
            endif 
          enddo; enddo
          if (counter > 0.0_pReal) FDot_inv = FDot_inv/counter
          do p=1_pInt,3_pInt; do o=1_pInt,3_pInt
            dFp_invdFdot(o,p,1:3,1:3) = Fpinv_rate(o,p)*FDot_inv
          enddo; enddo
          do p=1_pInt,3_pInt; do o=1_pInt,3_pInt
            dFedFdot(1:3,1:3,o,p) = math_transpose33(math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), &
                                    dFp_invdFdot(1:3,1:3,o,p)))
          enddo; enddo
          call constitutive_TandItsTangent(junk,dSdFe,crystallite_subFe0(1:3,1:3,g,i,e),g,i,e)            ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative
          dSdFdot = math_mul3333xx3333(dSdFe,dFedFdot)
          do p=1_pInt,3_pInt; do o=1_pInt,3_pInt
            crystallite_dPdF(1:3,1:3,o,p,g,i,e) = crystallite_dPdF(1:3,1:3,o,p,g,i,e) - &
                   (math_mul33x33(math_mul33x33(dFedFdot(1:3,1:3,o,p), &
                   math_Mandel6to33(crystallite_Tstar_v)),math_transpose33( &
                   crystallite_invFp(1:3,1:3,g,i,e))) + &                                                  ! dP/dFdot = dFe/dFdot * S * Fp^-T...
                   math_mul33x33(math_mul33x33(crystallite_subFe0(1:3,1:3,g,i,e), &
                   math_Mandel6to33(crystallite_Tstar_v)),math_transpose33(dFp_invdFdot(1:3,1:3,o,p))) &   !            + Fe * S * dFp^-T/dFdot...
                   + math_mul33x33(crystallite_subFe0(1:3,1:3,g,i,e), &
                   math_mul33x33(dSdFdot(1:3,1:3,o,p),math_transpose33(crystallite_invFp(1:3,1:3,g,i,e))))) !           + Fe * dS/dFdot * Fp^-T         
          enddo; enddo
    enddo; enddo; enddo 
  !$OMP END PARALLEL DO
  endif
endif                                                                                                   ! jacobian calculation
 
end subroutine crystallite_stressAndItsTangent


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state and Temperature with 4th order explicit Runge Kutta method 
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateRK4(gg,ii,ee)
 use prec, only:         pInt, &
                         pReal
 use numerics, only:     numerics_integrationMode
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_StateLoopDistribution
 use FEsolving, only:    FEsolving_execElem, & 
                         FEsolving_execIP
 use mesh, only:         mesh_element, &
                         mesh_NcpElems, &
                         mesh_maxNips
 use material, only:     homogenization_Ngrains, &
                         homogenization_maxNgrains
 use constitutive, only: constitutive_sizeDotState, &
                         constitutive_state, &
                         constitutive_subState0, &
                         constitutive_dotState, &
                         constitutive_RK4dotState, &
                         constitutive_collectDotState, &
                         constitutive_deltaState, &
                         constitutive_collectDeltaState, &
                         constitutive_dotTemperature, &
                         constitutive_microstructure
 
 implicit none
 
 real(pReal), dimension(4), parameter ::       timeStepFraction = [0.5_pReal, 0.5_pReal, 1.0_pReal, 1.0_pReal] ! factor giving the fraction of the original timestep used for Runge Kutta Integration
 real(pReal), dimension(4), parameter ::       weight = [1.0_pReal, 2.0_pReal, 2.0_pReal, 1.0_pReal]           ! weight of slope used for Runge Kutta integration
 
 !*** input variables ***!
 integer(pInt), optional, intent(in)::         ee, &                                                 ! element index
                                               ii, &                                                 ! integration point index
                                               gg                                                    ! grain index
 
 !*** output variables ***!
 
 !*** local variables ***!
 integer(pInt)                                 e, &                                                  ! element index in element loop
                                               i, &                                                  ! integration point index in ip loop
                                               g, &                                                  ! grain index in grain loop
                                               n, &
                                               mySizeDotState
 integer(pInt), dimension(2) ::                eIter                                                 ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                                              ! bounds for ip iteration
                                               gIter                                                 ! bounds for grain iteration
 real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                               RK4dotTemperature                                     ! evolution of Temperature of each grain for Runge Kutta integration
 logical                                       singleRun                                             ! flag indicating computation for single (g,i,e) triple
 
 
 if (present(ee) .and. present(ii) .and. present(gg)) then
   eIter = ee
   iIter(1:2,ee) = ii
   gIter(1:2,ee) = gg
   singleRun = .true.
 else
   eIter = FEsolving_execElem(1:2)
   do e = eIter(1),eIter(2)
     iIter(1:2,e) = FEsolving_execIP(1:2,e)
     gIter(1:2,e) = [ 1_pInt,homogenization_Ngrains(mesh_element(3,e))]
   enddo
   singleRun = .false.
 endif
 
 
 ! --- FIRST RUNGE KUTTA STEP ---
 
 RK4dotTemperature = 0.0_pReal                                                                             ! initialize Runge-Kutta dotTemperature
 !$OMP PARALLEL PRIVATE(mySizeDotState)
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     constitutive_RK4dotState(g,i,e)%p = 0.0_pReal                                                         ! initialize Runge-Kutta dotState
     if (crystallite_todo(g,i,e)) then
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                         crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
       crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                       crystallite_Temperature(g,i,e),g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
            .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
         if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         else                                                                                              ! if broken local...
           crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
         endif
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP END PARALLEL
 
 
 ! --- SECOND TO FOURTH RUNGE KUTTA STEP PLUS FINAL INTEGRATION ---
 
 do n = 1_pInt,4_pInt
 
   ! --- state update ---
 
   !$OMP PARALLEL PRIVATE(mySizeDotState)
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         if (n < 4) then
           constitutive_RK4dotState(g,i,e)%p = constitutive_RK4dotState(g,i,e)%p + weight(n)*constitutive_dotState(g,i,e)%p
           RK4dotTemperature(g,i,e) = RK4dotTemperature(g,i,e) + weight(n)*crystallite_dotTemperature(g,i,e)
         elseif (n == 4) then
           constitutive_dotState(g,i,e)%p = (constitutive_RK4dotState(g,i,e)%p + &
                                             weight(n)*constitutive_dotState(g,i,e)%p) / 6.0_pReal         ! use weighted RKdotState for final integration
           crystallite_dotTemperature(g,i,e) = (RK4dotTemperature(g,i,e) + weight(n)*crystallite_dotTemperature(g,i,e)) / 6.0_pReal
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO      
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                 + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e) * timeStepFraction(n)
         crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                 + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e) * timeStepFraction(n)
         if (n == 4) then                                                                                  ! final integration step
#ifndef _OPENMP
           if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
               .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                      .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
             mySizeDotState = constitutive_sizeDotState(g,i,e)
             write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
             write(6,*)
             write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
             write(6,*)
             write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
             write(6,*)
           endif
#endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO      
 
   
   ! --- state jump ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- update dependent states ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                           ! update dependent state variables to be consistent with basic states
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- stress integration ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e,timeStepFraction(n))                  ! fraction of original times step
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO      
 
   
   ! --- dot state and RK dot state---
 
   if (n < 4) then
     !$OMP DO
       do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
         if (crystallite_todo(g,i,e)) then
           call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                             crystallite_Temperature(g,i,e), timeStepFraction(n)*crystallite_subdt(g,i,e), & ! fraction of original timestep
                                             crystallite_subFrac, g,i,e)
           crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                           crystallite_Temperature(g,i,e),g,i,e)
         endif
       enddo; enddo; enddo
     !$OMP ENDDO
     !$OMP DO
       do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
         !$OMP FLUSH(crystallite_todo)
         if (crystallite_todo(g,i,e)) then
           if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                    ! NaN occured in dotState
                .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then         ! NaN occured in dotTemperature
             if (.not. crystallite_localPlasticity(g,i,e)) then                                            ! if broken non-local...
               !$OMP CRITICAL (checkTodo)
                 crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                     ! ...all non-locals skipped
               !$OMP END CRITICAL (checkTodo)
             else                                                                                          ! if broken local...
               crystallite_todo(g,i,e) = .false.                                                           ! ... skip this one next time
             endif
           endif
         endif
       enddo; enddo; enddo
     !$OMP ENDDO
   endif
   !$OMP END PARALLEL
   
 enddo
 
 
 ! --- SET CONVERGENCE FLAG ---
 
 do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
   if (crystallite_todo(g,i,e)) then
     crystallite_converged(g,i,e) = .true.                                                               ! if still "to do" then converged per definitionem
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
       !$OMP CRITICAL (distributionState)
         debug_StateLoopDistribution(4,numerics_integrationMode) = &
           debug_StateLoopDistribution(4,numerics_integrationMode) + 1_pInt
       !$OMP END CRITICAL (distributionState)
     endif
   endif
 enddo; enddo; enddo
 
 
 ! --- CHECK NONLOCAL CONVERGENCE ---
 
 if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
   if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) then                      ! any non-local not yet converged (or broken)...
     crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged
   endif
 endif
 
end subroutine crystallite_integrateStateRK4


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state and Temperature with 5th order Runge-Kutta Cash-Karp method with 
!> adaptive step size  (use 5th order solution to advance = "local extrapolation")
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateRKCK45(gg,ii,ee)
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_StateLoopDistribution
 use numerics, only:     rTol_crystalliteState, &
                         rTol_crystalliteTemperature, &
                         numerics_integrationMode
 use FEsolving, only:    FEsolving_execElem, & 
                         FEsolving_execIP
 use mesh, only:         mesh_element, &
                         mesh_NcpElems, &
                         mesh_maxNips
 use material, only:     homogenization_Ngrains, &
                         homogenization_maxNgrains
 use constitutive, only: constitutive_sizeDotState, &
                         constitutive_maxSizeDotState, &
                         constitutive_state, &
                         constitutive_aTolState, &
                         constitutive_subState0, &
                         constitutive_dotState, &
                         constitutive_RKCK45dotState, &
                         constitutive_collectDotState, &
                         constitutive_deltaState, &
                         constitutive_collectDeltaState, &
                         constitutive_dotTemperature, &
                         constitutive_microstructure
 
 implicit none
 !*** input variables ***!
 integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                               ii, &                     ! integration point index
                                               gg                        ! grain index
 !*** local variables ***!
 integer(pInt)                                 e, &                      ! element index in element loop
                                               i, &                      ! integration point index in ip loop
                                               g, &                      ! grain index in grain loop
                                               n, &                      ! stage index in integration stage loop
                                               mySizeDotState, &         ! size of dot State
                                               s                         ! state index
 integer(pInt), dimension(2) ::                eIter                     ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                  ! bounds for ip iteration
                                               gIter                     ! bounds for grain iteration
 real(pReal), dimension(6,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                               RKCK45dotTemperature      ! evolution of Temperature of each grain for Runge Kutta Cash Karp integration
 real(pReal), dimension(5,5) ::                a                         ! coefficients in Butcher tableau (used for preliminary integration in stages 2 to 6)
 real(pReal), dimension(6) ::                  b, db                     ! coefficients in Butcher tableau (used for final integration and error estimate)
 real(pReal), dimension(5) ::                  c                         ! coefficients in Butcher tableau (fractions of original time step in stages 2 to 6)
 real(pReal), dimension(constitutive_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                               stateResiduum, &          ! residuum from evolution in micrstructure
                                               relStateResiduum          ! relative residuum from evolution in microstructure
 real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                               temperatureResiduum, &    ! residuum from evolution in temperature
                                               relTemperatureResiduum    ! relative residuum from evolution in temperature
 logical                                       singleRun                 ! flag indicating computation for single (g,i,e) triple
 
 
 ! --- FILL BUTCHER TABLEAU ---
 
 a = 0.0_pReal
 b = 0.0_pReal
 db = 0.0_pReal
 c = 0.0_pReal
 
 a(1,1) = 0.2_pReal
 a(1,2) = 0.075_pReal
 a(2,2) = 0.225_pReal
 a(1,3) = 0.3_pReal
 a(2,3) = -0.9_pReal
 a(3,3) = 1.2_pReal
 a(1,4) = -11.0_pReal / 54.0_pReal
 a(2,4) = 2.5_pReal
 a(3,4) = -70.0_pReal / 27.0_pReal
 a(4,4) = 35.0_pReal / 27.0_pReal
 a(1,5) = 1631.0_pReal / 55296.0_pReal
 a(2,5) = 175.0_pReal / 512.0_pReal
 a(3,5) = 575.0_pReal / 13824.0_pReal
 a(4,5) = 44275.0_pReal / 110592.0_pReal
 a(5,5) = 253.0_pReal / 4096.0_pReal
 
 b(1) = 37.0_pReal / 378.0_pReal
 b(3) = 250.0_pReal / 621.0_pReal
 b(4) = 125.0_pReal / 594.0_pReal
 b(6) = 512.0_pReal / 1771.0_pReal
 
 db(1) = b(1) - 2825.0_pReal / 27648.0_pReal
 db(3) = b(3) - 18575.0_pReal / 48384.0_pReal
 db(4) = b(4) - 13525.0_pReal / 55296.0_pReal
 db(5) = - 277.0_pReal / 14336.0_pReal
 db(6) = b(6) - 0.25_pReal
 
 c(1) = 0.2_pReal
 c(2) = 0.3_pReal
 c(3) = 0.6_pReal
 c(4) = 1.0_pReal
 c(5) = 0.875_pReal
 
 
 ! --- LOOP ITERATOR FOR ELEMENT, GRAIN, IP ---
 
 if (present(ee) .and. present(ii) .and. present(gg)) then
   eIter = ee
   iIter(1:2,ee) = ii
   gIter(1:2,ee) = gg
   singleRun = .true.
 else
   eIter = FEsolving_execElem(1:2)
   do e = eIter(1),eIter(2)
     iIter(1:2,e) = FEsolving_execIP(1:2,e)
     gIter(1:2,e) = [1_pInt,homogenization_Ngrains(mesh_element(3,e))]
   enddo
   singleRun = .false.
 endif
 
 
 
 ! --- FIRST RUNGE KUTTA STEP ---
 
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(a,1x,i1)') '<< CRYST >> RUNGE KUTTA STEP',1
   !$OMP END CRITICAL (write2out)
 endif
 !$OMP PARALLEL PRIVATE(mySizeDotState)
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                         crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
       crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                       crystallite_Temperature(g,i,e),g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
            .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
         if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         else                                                                                              ! if broken local...
           crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
         endif
       endif
     endif
  enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP END PARALLEL
 
 
 ! --- SECOND TO SIXTH RUNGE KUTTA STEP ---
 
 do n = 1_pInt,5_pInt
 
   ! --- state update ---
   
   !$OMP PARALLEL PRIVATE(mySizeDotState)
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         constitutive_RKCK45dotState(n,g,i,e)%p = constitutive_dotState(g,i,e)%p                           ! store Runge-Kutta dotState
         RKCK45dotTemperature(n,g,i,e) = crystallite_dotTemperature(g,i,e)                                 ! store Runge-Kutta dotTemperature
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         if (n == 1) then                                                                                  ! NEED TO DO THE ADDITION IN THIS LENGTHY WAY BECAUSE OF PARALLELIZATION (CAN'T USE A REDUCTION CLAUSE ON A POINTER OR USER DEFINED TYPE)
           constitutive_dotState(g,i,e)%p = a(1,1) * constitutive_RKCK45dotState(1,g,i,e)%p
           crystallite_dotTemperature(g,i,e) = a(1,1) * RKCK45dotTemperature(1,g,i,e)            
         elseif (n == 2) then
           constitutive_dotState(g,i,e)%p = a(1,2) * constitutive_RKCK45dotState(1,g,i,e)%p &
                                          + a(2,2) * constitutive_RKCK45dotState(2,g,i,e)%p
           crystallite_dotTemperature(g,i,e) = a(1,2) * RKCK45dotTemperature(1,g,i,e) &
                                             + a(2,2) * RKCK45dotTemperature(2,g,i,e)
         elseif (n == 3) then
           constitutive_dotState(g,i,e)%p = a(1,3) * constitutive_RKCK45dotState(1,g,i,e)%p &
                                          + a(2,3) * constitutive_RKCK45dotState(2,g,i,e)%p &
                                          + a(3,3) * constitutive_RKCK45dotState(3,g,i,e)%p
           crystallite_dotTemperature(g,i,e) = a(1,3) * RKCK45dotTemperature(1,g,i,e) &
                                             + a(2,3) * RKCK45dotTemperature(2,g,i,e) &
                                             + a(3,3) * RKCK45dotTemperature(3,g,i,e)
         elseif (n == 4) then
           constitutive_dotState(g,i,e)%p = a(1,4) * constitutive_RKCK45dotState(1,g,i,e)%p &
                                          + a(2,4) * constitutive_RKCK45dotState(2,g,i,e)%p &
                                          + a(3,4) * constitutive_RKCK45dotState(3,g,i,e)%p &
                                          + a(4,4) * constitutive_RKCK45dotState(4,g,i,e)%p
           crystallite_dotTemperature(g,i,e) = a(1,4) * RKCK45dotTemperature(1,g,i,e) &
                                             + a(2,4) * RKCK45dotTemperature(2,g,i,e) &
                                             + a(3,4) * RKCK45dotTemperature(3,g,i,e) &
                                             + a(4,4) * RKCK45dotTemperature(4,g,i,e)
         elseif (n == 5) then
           constitutive_dotState(g,i,e)%p = a(1,5) * constitutive_RKCK45dotState(1,g,i,e)%p &
                                          + a(2,5) * constitutive_RKCK45dotState(2,g,i,e)%p &
                                          + a(3,5) * constitutive_RKCK45dotState(3,g,i,e)%p &
                                          + a(4,5) * constitutive_RKCK45dotState(4,g,i,e)%p &
                                          + a(5,5) * constitutive_RKCK45dotState(5,g,i,e)%p
           crystallite_dotTemperature(g,i,e) = a(1,5) * RKCK45dotTemperature(1,g,i,e) &
                                             + a(2,5) * RKCK45dotTemperature(2,g,i,e) &
                                             + a(3,5) * RKCK45dotTemperature(3,g,i,e) &
                                             + a(4,5) * RKCK45dotTemperature(4,g,i,e) &
                                             + a(5,5) * RKCK45dotTemperature(5,g,i,e)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                                       + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
         crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                        + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- state jump ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
   
   ! --- update dependent states ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                           ! update dependent state variables to be consistent with basic states
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- stress integration ---
   
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) =  crystallite_integrateStress(g,i,e,c(n))                                ! fraction of original time step
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO      
 
 
   ! --- dot state and RK dot state---
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a,1x,i1)') '<< CRYST >> Runge--Kutta step',n+1_pInt
   endif
#endif
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                           crystallite_Temperature(g,i,e), c(n)*crystallite_subdt(g,i,e), & ! fraction of original timestep
                                           crystallite_subFrac, g,i,e)
         crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                         crystallite_Temperature(g,i,e),g,i,e)
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
              .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
           if (.not. crystallite_localPlasticity(g,i,e)) then                                              ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                       ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                            ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP END PARALLEL
 
 enddo  
 
 
 ! --- STATE UPDATE WITH ERROR ESTIMATE FOR STATE AND TEMPERATURE ---
 
 relStateResiduum = 0.0_pReal
 relTemperatureResiduum = 0.0_pReal
 !$OMP PARALLEL PRIVATE(mySizeDotState)
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       mySizeDotState = constitutive_sizeDotState(g,i,e)
       constitutive_RKCK45dotState(6,g,i,e)%p = constitutive_dotState(g,i,e)%p                             ! store Runge-Kutta dotState
       RKCK45dotTemperature(6,g,i,e) = crystallite_dotTemperature(g,i,e)                                   ! store Runge-Kutta dotTemperature
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
       
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       mySizeDotState = constitutive_sizeDotState(g,i,e)
 
       ! --- absolute residuum in state and temperature ---
       ! NEED TO DO THE ADDITION IN THIS LENGTHY WAY BECAUSE OF PARALLELIZATION 
       ! CAN'T USE A REDUCTION CLAUSE ON A POINTER OR USER DEFINED TYPE
       
       stateResiduum(1:mySizeDotState,g,i,e) = &
                 (   db(1) * constitutive_RKCK45dotState(1,g,i,e)%p(1:mySizeDotState)  &
                   + db(2) * constitutive_RKCK45dotState(2,g,i,e)%p(1:mySizeDotState)  &
                   + db(3) * constitutive_RKCK45dotState(3,g,i,e)%p(1:mySizeDotState)  &
                   + db(4) * constitutive_RKCK45dotState(4,g,i,e)%p(1:mySizeDotState)  &
                   + db(5) * constitutive_RKCK45dotState(5,g,i,e)%p(1:mySizeDotState)  &
                   + db(6) * constitutive_RKCK45dotState(6,g,i,e)%p(1:mySizeDotState)) &
                                    * crystallite_subdt(g,i,e)
       temperatureResiduum(g,i,e) = (  db(1) * RKCK45dotTemperature(1,g,i,e) &
                                     + db(2) * RKCK45dotTemperature(2,g,i,e) &
                                     + db(3) * RKCK45dotTemperature(3,g,i,e) &
                                     + db(4) * RKCK45dotTemperature(4,g,i,e) &
                                     + db(5) * RKCK45dotTemperature(5,g,i,e) &
                                     + db(6) * RKCK45dotTemperature(6,g,i,e)) &
                                    * crystallite_subdt(g,i,e)
 
       ! --- dot state and dot temperature ---
 
       constitutive_dotState(g,i,e)%p = b(1) * constitutive_RKCK45dotState(1,g,i,e)%p &
                                      + b(2) * constitutive_RKCK45dotState(2,g,i,e)%p &
                                      + b(3) * constitutive_RKCK45dotState(3,g,i,e)%p &
                                      + b(4) * constitutive_RKCK45dotState(4,g,i,e)%p &
                                      + b(5) * constitutive_RKCK45dotState(5,g,i,e)%p &
                                      + b(6) * constitutive_RKCK45dotState(6,g,i,e)%p
       crystallite_dotTemperature(g,i,e) = b(1) * RKCK45dotTemperature(1,g,i,e) &
                                         + b(2) * RKCK45dotTemperature(2,g,i,e) &
                                         + b(3) * RKCK45dotTemperature(3,g,i,e) &
                                         + b(4) * RKCK45dotTemperature(4,g,i,e) &
                                         + b(5) * RKCK45dotTemperature(5,g,i,e) &
                                         + b(6) * RKCK45dotTemperature(6,g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 ! --- state and temperature update ---      
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       mySizeDotState = constitutive_sizeDotState(g,i,e)
       constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                                     + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
       crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                      + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 ! --- relative residui and state convergence ---      
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       mySizeDotState = constitutive_sizeDotState(g,i,e)
       forall (s = 1_pInt:mySizeDotState, abs(constitutive_state(g,i,e)%p(s)) > 0.0_pReal) &
         relStateResiduum(s,g,i,e) = stateResiduum(s,g,i,e) / constitutive_state(g,i,e)%p(s)
       if (crystallite_Temperature(g,i,e) > 0) then
         relTemperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) / crystallite_Temperature(g,i,e)
       endif
       !$OMP FLUSH(relStateResiduum,relTemperatureResiduum)
       
       crystallite_todo(g,i,e) = &
           ( all(      abs(relStateResiduum(:,g,i,e)) < rTol_crystalliteState &
                  .or. abs(stateResiduum(1:mySizeDotState,g,i,e)) < constitutive_aTolState(g,i,e)%p(1:mySizeDotState) ) &
            .and. abs(relTemperatureResiduum(g,i,e)) < rTol_crystalliteTemperature )
            
#ifndef _OPENMP
       if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt&
           .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                  .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
         write(6,'(a,i8,1x,i3,1x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
         write(6,*)
         write(6,'(a,/,(12x,12(f12.1,1x)))') '<< CRYST >> absolute residuum tolerance', &
                                         stateResiduum(1:mySizeDotState,g,i,e) / constitutive_aTolState(g,i,e)%p(1:mySizeDotState)
         write(6,*)
         write(6,'(a,/,(12x,12(f12.1,1x)))') '<< CRYST >> relative residuum tolerance', &
                                               relStateResiduum(1:mySizeDotState,g,i,e) / rTol_crystalliteState
         write(6,*)
         write(6,'(a,/,(12x,12(e12.5,1x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
         write(6,*)
         write(6,'(a,/,(12x,12(e12.5,1x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
         write(6,*)
       endif
#endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
       
 
 ! --- STATE JUMP ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 
 ! --- UPDATE DEPENDENT STATES IF RESIDUUM BELOW TOLERANCE ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                        crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                             ! update dependent state variables to be consistent with basic states
     endif
  enddo; enddo; enddo
 !$OMP ENDDO
 
 
 ! --- FINAL STRESS INTEGRATION STEP IF RESIDUUM BELOW TOLERANCE ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then              ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 
 ! --- SET CONVERGENCE FLAG ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       crystallite_converged(g,i,e) = .true.                                                               ! if still "to do" then converged per definitionem
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (distributionState)
           debug_StateLoopDistribution(6,numerics_integrationMode) = &
             debug_StateLoopDistribution(6,numerics_integrationMode) + 1_pInt
         !$OMP END CRITICAL (distributionState)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 !$OMP END PARALLEL
 
 
 ! --- nonlocal convergence check ---
 
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), ' grains converged'
   write(6,*)
   !$OMP END CRITICAL (write2out)
 endif
 if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
   if ( any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) then                     ! any non-local not yet converged (or broken)...
     crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged
   endif
   
 endif
 
end subroutine crystallite_integrateStateRKCK45


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state and Temperature with 1st order Euler method with adaptive step size
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateAdaptiveEuler(gg,ii,ee)
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_StateLoopDistribution
 use numerics, only:     rTol_crystalliteState, &
                         rTol_crystalliteTemperature, &
                         numerics_integrationMode
 use FEsolving, only:    FEsolving_execElem, & 
                         FEsolving_execIP
 use mesh, only:         mesh_element, &
                         mesh_NcpElems, &
                         mesh_maxNips
 use material, only:     homogenization_Ngrains, &
                         homogenization_maxNgrains
 use constitutive, only: constitutive_sizeDotState, &
                         constitutive_maxSizeDotState, &
                         constitutive_state, &
                         constitutive_aTolState, &
                         constitutive_subState0, &
                         constitutive_dotState, &
                         constitutive_collectDotState, &
                         constitutive_dotTemperature, &
                         constitutive_microstructure
 
 implicit none
 
 !*** input variables ***!
 integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                               ii, &                     ! integration point index
                                               gg                        ! grain index
 !*** local variables ***!
 integer(pInt)                                 e, &                      ! element index in element loop
                                               i, &                      ! integration point index in ip loop
                                               g, &                      ! grain index in grain loop
                                               mySizeDotState, &         ! size of dot State
                                               s                         ! state index
 integer(pInt), dimension(2) ::                eIter                     ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                  ! bounds for ip iteration
                                               gIter                     ! bounds for grain iteration
 real(pReal), dimension(constitutive_maxSizeDotState,homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                               stateResiduum, &          ! residuum from evolution in micrstructure
                                               relStateResiduum          ! relative residuum from evolution in microstructure
 real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                               temperatureResiduum, &    ! residuum from evolution in temperature
                                               relTemperatureResiduum    ! relative residuum from evolution in temperature
 logical                                       singleRun                 ! flag indicating computation for single (g,i,e) triple
 
 
 ! --- LOOP ITERATOR FOR ELEMENT, GRAIN, IP ---
 
 if (present(ee) .and. present(ii) .and. present(gg)) then
   eIter = ee
   iIter(1:2,ee) = ii
   gIter(1:2,ee) = gg
   singleRun = .true.
 else
   eIter = FEsolving_execElem(1:2)
   do e = eIter(1),eIter(2)
     iIter(1:2,e) = FEsolving_execIP(1:2,e)
     gIter(1:2,e) = [1_pInt,homogenization_Ngrains(mesh_element(3,e))]
   enddo
   singleRun = .false.
 endif
 
 
 stateResiduum = 0.0_pReal
 
 !$OMP PARALLEL PRIVATE(mySizeDotState)
 
 if (numerics_integrationMode == 1_pInt) then
 
   ! --- DOT STATE AND TEMPERATURE (EULER INTEGRATION) ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then  
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                           crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
         crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                         crystallite_Temperature(g,i,e),g,i,e)
       endif
    enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then  
         if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
              .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                              ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
           endif
         endif
       endif
    enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- STATE UPDATE (EULER INTEGRATION) ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         stateResiduum(1:mySizeDotState,g,i,e) = - 0.5_pReal * constitutive_dotState(g,i,e)%p * crystallite_subdt(g,i,e)  ! contribution to absolute residuum in state and temperature
         temperatureResiduum(g,i,e) = - 0.5_pReal * crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
         constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_state(g,i,e)%p(1:mySizeDotState) &
                                                       + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
         crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                        + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
       endif
     enddo; enddo; enddo
   !$OMP ENDDO      
 
 
   ! --- STATE JUMP ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- UPDATE DEPENDENT STATES (EULER INTEGRATION) ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                             ! update dependent state variables to be consistent with basic states
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 endif
 
 ! --- STRESS INTEGRATION (EULER INTEGRATION) ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! if broken non-local...
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO      
 
 
 if (numerics_integrationMode == 1_pInt) then
 
   ! --- DOT STATE AND TEMPERATURE (HEUN METHOD) ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                           crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
         crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                         crystallite_Temperature(g,i,e),g,i,e)
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e)) then
         if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                      ! NaN occured in dotState
              .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then           ! NaN occured in dotTemperature
           if (.not. crystallite_localPlasticity(g,i,e)) then                                              ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                       ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                            ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
           endif
         endif      
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- ERROR ESTIMATE FOR STATE AND TEMPERATURE (HEUN METHOD) ---
 
   !$OMP SINGLE
   relStateResiduum = 0.0_pReal
   relTemperatureResiduum = 0.0_pReal
   !$OMP END SINGLE
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         
         
         ! --- contribution of heun step to absolute residui ---
         
         stateResiduum(1:mySizeDotState,g,i,e) = stateResiduum(1:mySizeDotState,g,i,e) &
                                               + 0.5_pReal * constitutive_dotState(g,i,e)%p * crystallite_subdt(g,i,e) ! contribution to absolute residuum in state and temperature
         temperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) &
                                    + 0.5_pReal * crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
         !$OMP FLUSH(stateResiduum,temperatureResiduum)
 
         ! --- relative residui ---      
 
         forall (s = 1_pInt:mySizeDotState, abs(constitutive_state(g,i,e)%p(s)) > 0.0_pReal) &
           relStateResiduum(s,g,i,e) = stateResiduum(s,g,i,e) / constitutive_state(g,i,e)%p(s)
         if (crystallite_Temperature(g,i,e) > 0_pInt) then
           relTemperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) / crystallite_Temperature(g,i,e)
         endif
         !$OMP FLUSH(relStateResiduum,relTemperatureResiduum)
 
#ifndef _OPENMP        
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
           write(6,*)
           write(6,'(a,/,(12x,12(f12.1,1x)))') '<< CRYST >> absolute residuum tolerance', &
                                           stateResiduum(1:mySizeDotState,g,i,e) / constitutive_aTolState(g,i,e)%p(1:mySizeDotState)
           write(6,*)
           write(6,'(a,/,(12x,12(f12.1,1x)))') '<< CRYST >> relative residuum tolerance', &
                                                 relStateResiduum(1:mySizeDotState,g,i,e) / rTol_crystalliteState
           write(6,*)
           write(6,'(a,/,(12x,12(e12.5,1x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState) &
                                                    - 2.0_pReal * stateResiduum(1:mySizeDotState,g,i,e) / crystallite_subdt(g,i,e)  ! calculate former dotstate from higher order solution and state residuum
           write(6,*)
           write(6,'(a,/,(12x,12(e12.5,1x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
           write(6,*)
         endif
#endif
         
         ! --- converged ? ---
 
         if ( all(     abs(relStateResiduum(:,g,i,e)) < rTol_crystalliteState &
                  .or. abs(stateResiduum(1:mySizeDotState,g,i,e)) < constitutive_aTolState(g,i,e)%p(1:mySizeDotState)) &
             .and. abs(relTemperatureResiduum(g,i,e)) < rTol_crystalliteTemperature ) then        
           crystallite_converged(g,i,e) = .true.                                                             ! ... converged per definitionem
           if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionState)
               debug_StateLoopDistribution(2,numerics_integrationMode) = &
                 debug_StateLoopDistribution(2,numerics_integrationMode) + 1_pInt
             !$OMP END CRITICAL (distributionState)
           endif
         endif
 
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 elseif (numerics_integrationMode > 1) then ! stiffness calculation
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e)) then
         crystallite_converged(g,i,e) = .true.                                                               ! ... converged per definitionem
         if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
           !$OMP CRITICAL (distributionState)
             debug_StateLoopDistribution(2,numerics_integrationMode) = &
               debug_StateLoopDistribution(2,numerics_integrationMode) + 1_pInt
           !$OMP END CRITICAL (distributionState)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 endif
 
 !$OMP END PARALLEL
 
 
 ! --- NONLOCAL CONVERGENCE CHECK ---
 
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
   !$OMP CRITICAL (write2out)
   write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), ' grains converged'
   write(6,*)
   !$OMP END CRITICAL (write2out)
 endif
 if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
   if ( any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) then                     ! any non-local not yet converged (or broken)...
     crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged
   endif
 endif
 
end subroutine crystallite_integrateStateAdaptiveEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state and Temperature with 1st order explicit Euler method
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateEuler(gg,ii,ee)
 use numerics, only:     numerics_integrationMode, &
                         numerics_timeSyncing
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_StateLoopDistribution
 use FEsolving, only:    FEsolving_execElem, & 
                         FEsolving_execIP
 use mesh, only:         mesh_element, &
                         mesh_NcpElems
 use material, only:     homogenization_Ngrains
 use constitutive, only: constitutive_sizeDotState, &
                         constitutive_state, &
                         constitutive_subState0, &
                         constitutive_dotState, &
                         constitutive_collectDotState, &
                         constitutive_dotTemperature, &
                         constitutive_microstructure
 
 implicit none
 !*** input variables ***!
 integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                               ii, &                     ! integration point index
                                               gg                        ! grain index
 !*** local variables ***!
 integer(pInt)                                 e, &                      ! element index in element loop
                                               i, &                      ! integration point index in ip loop
                                               g, &                      ! grain index in grain loop
                                               mySizeDotState
 integer(pInt), dimension(2) ::                eIter                     ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                  ! bounds for ip iteration
                                               gIter                     ! bounds for grain iteration
 logical                                       singleRun                 ! flag indicating computation for single (g,i,e) triple
 
 
 if (present(ee) .and. present(ii) .and. present(gg)) then
   eIter = ee
   iIter(1:2,ee) = ii
   gIter(1:2,ee) = gg
   singleRun = .true.
 else
   eIter = FEsolving_execElem(1:2)
   do e = eIter(1),eIter(2)
     iIter(1:2,e) = FEsolving_execIP(1:2,e)
     gIter(1:2,e) = [1_pInt,homogenization_Ngrains(mesh_element(3,e))]
   enddo
   singleRun = .false.
 endif
 
 
 !$OMP PARALLEL
 
 if (numerics_integrationMode == 1_pInt) then
 
   ! --- DOT STATE AND TEMPERATURE ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                           crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
         crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                         crystallite_Temperature(g,i,e),g,i,e)
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
              .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
           if (.not. crystallite_localPlasticity(g,i,e) .and. .not. numerics_timeSyncing) then              ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           else                                                                                              ! if broken local...
             crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- UPDATE STATE AND TEMPERATURE ---
 
   !$OMP DO PRIVATE(mySizeDotState)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_state(g,i,e)%p(1:mySizeDotState) &
                                                       + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
         crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                        + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                     .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> update state at el ip g ',e,i,g
           write(6,*)
           write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
           write(6,*)
           write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
           write(6,*)
         endif
#endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- STATE JUMP ---
   
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e) &                ! if broken non-local...
             .and. .not. numerics_timeSyncing) then
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- UPDATE DEPENDENT STATES ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                            ! update dependent state variables to be consistent with basic states
       endif
    enddo; enddo; enddo
   !$OMP ENDDO
 
 endif
 
 
 ! --- STRESS INTEGRATION ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
       !$OMP FLUSH(crystallite_todo)
       if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e) &                  ! if broken non-local...
           .and. .not. numerics_timeSyncing) then
         !$OMP CRITICAL (checkTodo)
           crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                           ! ...all non-locals skipped
         !$OMP END CRITICAL (checkTodo)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 
 ! --- SET CONVERGENCE FLAG ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
       crystallite_converged(g,i,e) = .true.                                                               ! if still "to do" then converged per definitionem
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         !$OMP CRITICAL (distributionState)
           debug_StateLoopDistribution(1,numerics_integrationMode) = &
             debug_StateLoopDistribution(1,numerics_integrationMode) + 1_pInt
         !$OMP END CRITICAL (distributionState)
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 !$OMP END PARALLEL
 
 
 ! --- CHECK NON-LOCAL CONVERGENCE ---
 
 if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
   if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity) &                          ! any non-local not yet converged (or broken)...
       .and. .not. numerics_timeSyncing) then
     crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                       ! ...restart all non-local as not converged
   endif
 endif
 
end subroutine crystallite_integrateStateEuler


!--------------------------------------------------------------------------------------------------
!> @brief integrate stress, state and Temperature with adaptive 1st order explicit Euler method  
!> using Fixed Point Iteration to adapt the stepsize  
!--------------------------------------------------------------------------------------------------
subroutine crystallite_integrateStateFPI(gg,ii,ee)
 use debug, only:        debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_level,&
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_StateLoopDistribution
 use numerics, only:     nState, &
                         numerics_integrationMode, &
                         rTol_crystalliteState, &
                         rTol_crystalliteTemperature
 use FEsolving, only:    FEsolving_execElem, & 
                         FEsolving_execIP
 use mesh, only:         mesh_element, &
                         mesh_NcpElems
 use material, only:     homogenization_Ngrains
 use constitutive, only: constitutive_subState0, &
                         constitutive_state, &
                         constitutive_sizeDotState, &
                         constitutive_maxSizeDotState, &
                         constitutive_dotState, &
                         constitutive_collectDotState, &
                         constitutive_dotTemperature, &
                         constitutive_microstructure, &
                         constitutive_previousDotState, &
                         constitutive_previousDotState2, &
                         constitutive_aTolState
 
 implicit none
 integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                               ii, &                     ! integration point index
                                               gg                        ! grain index
 
 !*** local variables ***!
 integer(pInt)                                 NiterationState, &        ! number of iterations in state loop
                                               e, &                      ! element index in element loop
                                               i, &                      ! integration point index in ip loop
                                               g, &                      ! grain index in grain loop
                                               mySizeDotState
 integer(pInt), dimension(2) ::                eIter                     ! bounds for element iteration
 integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                  ! bounds for ip iteration
                                               gIter                     ! bounds for grain iteration
 real(pReal)                                   dot_prod12, &
                                               dot_prod22, &
                                               stateDamper, &            ! damper for integration of state
                                               temperatureResiduum
 real(pReal), dimension(constitutive_maxSizeDotState) :: &
                                               stateResiduum, &
                                               tempState
 logical                                       singleRun                 ! flag indicating computation for single (g,i,e) triple
 
 singleRun = present(ee) .and. present(ii) .and. present(gg)
 if (singleRun) then
   eIter = ee
   iIter(1:2,ee) = ii
   gIter(1:2,ee) = gg
 else
   eIter = FEsolving_execElem(1:2)
   do e = eIter(1),eIter(2)
     iIter(1:2,e) = FEsolving_execIP(1:2,e)
     gIter(1:2,e) = [1_pInt,homogenization_Ngrains(mesh_element(3,e))]
   enddo
 endif
 
 
 ! --+>> PREGUESS FOR STATE <<+--
 
 !$OMP PARALLEL
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     constitutive_previousDotState(g,i,e)%p = 0.0_pReal
     constitutive_previousDotState2(g,i,e)%p = 0.0_pReal
  enddo; enddo; enddo
 !$OMP ENDDO
 
 
 ! --- DOT STATES ---
 
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                         crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
       crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                       crystallite_Temperature(g,i,e),g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 !$OMP DO
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     !$OMP FLUSH(crystallite_todo)
     if (crystallite_todo(g,i,e)) then
       if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
            .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
         if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken is a non-local...
           !$OMP CRITICAL (checkTodo)
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals done (and broken)
           !$OMP END CRITICAL (checkTodo)
         else                                                                                              ! broken one was local...
           crystallite_todo(g,i,e) = .false.                                                               ! ... done (and broken)
         endif
       endif
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 
 ! --- UPDATE STATE AND TEMPERATURE ---
 
 !$OMP DO PRIVATE(mySizeDotState)
   do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
     if (crystallite_todo(g,i,e)) then
       mySizeDotState = constitutive_sizeDotState(g,i,e)
       constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                                     + constitutive_dotState(g,i,e)%p * crystallite_subdt(g,i,e)
       crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                      + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
     endif
   enddo; enddo; enddo
 !$OMP ENDDO
 
 !$OMP END PARALLEL
 
 
 ! --+>> STATE LOOP <<+--
 
 NiterationState = 0_pInt
 do while (any(crystallite_todo .and. .not. crystallite_converged) .and. NiterationState < nState )          ! convergence loop for crystallite
   NiterationState = NiterationState + 1_pInt
   
   !$OMP PARALLEL
 
   ! --- UPDATE DEPENDENT STATES ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe(1:3,1:3,g,i,e), &
                                          crystallite_Fp(1:3,1:3,g,i,e), g, i, e)                            ! update dependent state variables to be consistent with basic states
       endif
       constitutive_previousDotState2(g,i,e)%p = constitutive_previousDotState(g,i,e)%p                      ! remember previous dotState
       constitutive_previousDotState(g,i,e)%p = constitutive_dotState(g,i,e)%p                               ! remember current  dotState
     enddo; enddo; enddo
   !$OMP ENDDO
   
   
   ! --- STRESS INTEGRATION ---
   
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         crystallite_todo(g,i,e) = crystallite_integrateStress(g,i,e)
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e) .and. .not. crystallite_localPlasticity(g,i,e)) then            ! broken non-local... 
           !$OMP CRITICAL (checkTodo) 
             crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ... then all non-locals skipped
           !$OMP END CRITICAL (checkTodo)
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   !$OMP SINGLE
   !$OMP CRITICAL (write2out)
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
     write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo after stress integration'
   endif
   !$OMP END CRITICAL (write2out)
   !$OMP END SINGLE
 
 
   ! --- DOT STATE AND TEMPERATURE ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                           crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_subFrac, g,i,e)
         crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                         crystallite_Temperature(g,i,e),g,i,e)
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
         if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
              .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
           crystallite_todo(g,i,e) = .false.                                                                 ! ... skip me next time
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if me is non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- UPDATE STATE AND TEMPERATURE ---
 
   !$OMP DO PRIVATE(dot_prod12,dot_prod22,statedamper,mySizeDotState,stateResiduum,temperatureResiduum,tempState)
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       if (crystallite_todo(g,i,e) .and. .not. crystallite_converged(g,i,e)) then
 
         ! --- state damper ---
         
         dot_prod12 = dot_product( constitutive_dotState(g,i,e)%p - constitutive_previousDotState(g,i,e)%p, &
                           constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p )
         dot_prod22 = dot_product( constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p, &
                                   constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p )
         if (      dot_prod22 > 0.0_pReal &
             .and. (     dot_prod12 < 0.0_pReal &
                    .or. dot_product(constitutive_dotState(g,i,e)%p, constitutive_previousDotState(g,i,e)%p) < 0.0_pReal) ) then
           statedamper = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
         else
           statedamper = 1.0_pReal
         endif
         
         ! --- get residui ---
         
         mySizeDotState = constitutive_sizeDotState(g,i,e)
         stateResiduum(1:mySizeDotState) = constitutive_state(g,i,e)%p(1:mySizeDotState) &
                                 - constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                 - (constitutive_dotState(g,i,e)%p * statedamper &
                                    + constitutive_previousDotState(g,i,e)%p * (1.0_pReal - statedamper)) * crystallite_subdt(g,i,e)
         temperatureResiduum = crystallite_Temperature(g,i,e) &
                             - crystallite_subTemperature0(g,i,e) &
                             - crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
         
         ! --- correct state and temperature with residuum ---
 
         tempState(1:mySizeDotState) = constitutive_state(g,i,e)%p(1:mySizeDotState) - stateResiduum(1:mySizeDotState) ! need to copy to local variable, since we cant flush a pointer in openmp
         crystallite_Temperature(g,i,e) = crystallite_Temperature(g,i,e) - temperatureResiduum
         !$OMP FLUSH(crystallite_Temperature)
#ifndef _OPENMP
         if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> update state at el ip g ',e,i,g
           write(6,*)
           write(6,'(a,f6.1)') '<< CRYST >> statedamper ',statedamper
           write(6,*)
           write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> state residuum',stateResiduum(1:mySizeDotState)
           write(6,*)
           write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> new state',tempState(1:mySizeDotState)
           write(6,*)
         endif
#endif
 
         ! --- store corrected dotState --- (cannot do this before state update, because not sure how to flush pointers in openmp)
 
         constitutive_dotState(g,i,e)%p = constitutive_dotState(g,i,e)%p * statedamper &
                                        + constitutive_previousDotState(g,i,e)%p * (1.0_pReal - statedamper)
 
 
         ! --- converged ? ---
 
         if ( all(    abs(stateResiduum(1:mySizeDotState)) < constitutive_aTolState(g,i,e)%p(1:mySizeDotState) &
                 .or. abs(stateResiduum(1:mySizeDotState)) < rTol_crystalliteState &
                                                           * abs(tempState(1:mySizeDotState)) ) &
             .and. (abs(temperatureResiduum) < rTol_crystalliteTemperature * crystallite_Temperature(g,i,e) &
                    .or. crystallite_Temperature(g,i,e) == 0.0_pReal) ) then        
           crystallite_converged(g,i,e) = .true.                                                             ! ... converged per definitionem
           if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
             !$OMP CRITICAL (distributionState)
               debug_StateLoopDistribution(NiterationState,numerics_integrationMode) = &
                 debug_StateLoopDistribution(NiterationState,numerics_integrationMode) + 1_pInt
             !$OMP END CRITICAL (distributionState)
           endif
         endif
         constitutive_state(g,i,e)%p(1:mySizeDotState) = tempState(1:mySizeDotState)  ! copy local backup to global pointer
         
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
 
   ! --- STATE JUMP ---
 
   !$OMP DO
     do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
       !$OMP FLUSH(crystallite_todo)
       if (crystallite_todo(g,i,e) .and. crystallite_converged(g,i,e)) then                                  ! converged and still alive...
         crystallite_todo(g,i,e) = crystallite_stateJump(g,i,e)                                         
         !$OMP FLUSH(crystallite_todo)
         if (.not. crystallite_todo(g,i,e)) then                                                             ! if state jump fails, then convergence is broken
           crystallite_converged(g,i,e) = .false.
           if (.not. crystallite_localPlasticity(g,i,e)) then                                                ! if broken non-local...
             !$OMP CRITICAL (checkTodo)
               crystallite_todo = crystallite_todo .and. crystallite_localPlasticity                         ! ...all non-locals skipped
             !$OMP END CRITICAL (checkTodo)
           endif
         endif
       endif
     enddo; enddo; enddo
   !$OMP ENDDO
 
   !$OMP END PARALLEL
 
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL(write2out)
     write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), &
                               ' grains converged after state integration no. ', NiterationState
     write(6,*)
     !$OMP END CRITICAL(write2out)
   endif
 
   
   ! --- NON-LOCAL CONVERGENCE CHECK ---
 
   if (.not. singleRun) then                                                                               ! if not requesting Integration of just a single IP   
     if (any(.not. crystallite_converged .and. .not. crystallite_localPlasticity)) then                    ! any non-local not yet converged (or broken)...
       crystallite_converged = crystallite_converged .and. crystallite_localPlasticity                     ! ...restart all non-local as not converged
     endif
   endif
   
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt) then
     !$OMP CRITICAL(write2out)
     write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_converged(:,:,:)),' grains converged after non-local check'
     write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo after state integration no. ',&
                                 NiterationState
     write(6,*)
     !$OMP END CRITICAL(write2out)
   endif
   
 enddo                                                                                           ! crystallite convergence loop  
 
end subroutine crystallite_integrateStateFPI


!--------------------------------------------------------------------------------------------------
!> @brief calculates a jump in the state according to the current state and the current stress
!--------------------------------------------------------------------------------------------------
function crystallite_stateJump(g,i,e)
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g
 use FEsolving, only:    FEsolving_execElem, & 
                         FEsolving_execIP
 use mesh, only:         mesh_element, &
                         mesh_NcpElems
 use material, only:     homogenization_Ngrains
 use constitutive, only: constitutive_sizeDotState, &
                         constitutive_state, &
                         constitutive_deltaState, &
                         constitutive_collectDeltaState, &
                         constitutive_microstructure
 
 implicit none
 integer(pInt), intent(in):: e, &                      ! element index
                             i, &                      ! integration point index
                             g                         ! grain index
 
 !*** output variables ***!
 logical                     crystallite_stateJump
 
 !*** local variables ***!
 integer(pInt)               mySizeDotState
 
 
 crystallite_stateJump = .false.
 
 call constitutive_collectDeltaState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Temperature(g,i,e), g,i,e)
 
 mySizeDotState = constitutive_sizeDotState(g,i,e)
 if (any(constitutive_deltaState(g,i,e)%p(1:mySizeDotState) /= constitutive_deltaState(g,i,e)%p(1:mySizeDotState))) then
   return
 endif
 
 constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_state(g,i,e)%p(1:mySizeDotState) &
                                               + constitutive_deltaState(g,i,e)%p(1:mySizeDotState)
 
#ifndef _OPENMP
 if (any(constitutive_deltaState(g,i,e)%p(1:mySizeDotState) /= 0.0_pReal) &
     .and. iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
     .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> update state at el ip g ',e,i,g
   write(6,*)
   write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> deltaState', constitutive_deltaState(g,i,e)%p(1:mySizeDotState)
   write(6,*)
   write(6,'(a,/,(12x,12(e12.6,1x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
   write(6,*)
 endif
#endif
 
 crystallite_stateJump = .true.
 
end function crystallite_stateJump



!--------------------------------------------------------------------------------------------------
!> @brief calculation of stress (P) with time integration based on a residuum in Lp and 
!> intermediate acceleration of the Newton-Raphson correction  
!--------------------------------------------------------------------------------------------------
function crystallite_integrateStress(&
      g,&          ! grain number
      i,&          ! integration point number
      e,&          ! element number
      timeFraction &
      )
 use prec, only:         pLongInt
 use numerics, only:     nStress, &
                         aTol_crystalliteStress, &
                         rTol_crystalliteStress, &
                         iJacoLpresiduum, &
                         numerics_integrationMode
 use debug, only:        debug_level, &
                         debug_crystallite, &
                         debug_levelBasic, &
                         debug_levelExtensive, &
                         debug_levelSelective, &
                         debug_e, &
                         debug_i, &
                         debug_g, &
                         debug_cumLpCalls, &
                         debug_cumLpTicks, &
                         debug_StressLoopDistribution
 use constitutive, only: constitutive_LpAndItsTangent, &
                         constitutive_TandItsTangent, &
                         constitutive_homogenizedC
 use math, only:         math_mul33x33, &
                         math_mul33xx33, &
                         math_mul66x6, &
                         math_mul99x99, &
                         math_transpose33, &
                         math_inv33, &
                         math_invert33, &
                         math_invert, &
                         math_det33, &
                         math_norm33, &
                         math_I3, &
                         math_identity2nd, &
                         math_Mandel66to3333, &
                         math_Mandel6to33, &
                         math_Mandel33to6, &
                         math_Plain3333to99, &
                         math_Plain33to9, &
                         math_Plain9to33
 
 implicit none
 integer(pInt), intent(in)::         e, &                          ! element index
                                     i, &                          ! integration point index
                                     g                             ! grain index
 real(pReal), optional, intent(in) :: timeFraction                 ! fraction of timestep
 
 !*** output variables ***!
 logical                             crystallite_integrateStress   ! flag indicating if integration suceeded
 
 !*** local variables ***!
 real(pReal), dimension(3,3)::       Fg_new, &                                                       ! deformation gradient at end of timestep
                                     Fp_current, &                                                   ! plastic deformation gradient at start of timestep
                                     Fp_new, &                                                       ! plastic deformation gradient at end of timestep
                                     Fe_new, &                                                       ! elastic deformation gradient at end of timestep
                                     invFp_new, &                                                    ! inverse of Fp_new
                                     invFp_current, &                                                ! inverse of Fp_current
                                     Lpguess, &                                                      ! current guess for plastic velocity gradient
                                     Lpguess_old, &                                                  ! known last good guess for plastic velocity gradient
                                     Lp_constitutive, &                                              ! plastic velocity gradient resulting from constitutive law
                                     residuum, &                                                     ! current residuum of plastic velocity gradient
                                     residuum_old, &                                                 ! last residuum of plastic velocity gradient
                                     deltaLp, &                                                      ! direction of next guess
                                     Tstar,&                                                         ! 2nd Piola-Kirchhoff Stress
                                     A,&
                                     B, &
                                     Fe                                                              ! elastic deformation gradient
 real(pReal), dimension(6)::         Tstar_v                                                         ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
 real(pReal), dimension(9)::         work                                                            ! needed for matrix inversion by LAPACK
 integer(pInt), dimension(9) ::      ipiv                                                            ! needed for matrix inversion by LAPACK
 real(pReal), dimension(9,9) ::      dLp_dT_constitutive, &                                          ! partial derivative of plastic velocity gradient calculated by constitutive law
                                     dT_dFe_constitutive, &                                          ! partial derivative of 2nd Piola-Kirchhoff stress calculated by constitutive law
                                     dFe_dLp, &                                                      ! partial derivative of elastic deformation gradient
                                     dR_dLp, &                                                       ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                     dR_dLp2                                                         ! working copy of dRdLp
 real(pReal), dimension(3,3,3,3)::   dT_dFe3333, &                                                   ! partial derivative of 2nd Piola-Kirchhoff stress
                                     dFe_dLp3333                                                     ! partial derivative of elastic deformation gradient
 real(pReal)                         p_hydro, &                                                      ! volumetric part of 2nd Piola-Kirchhoff Stress
                                     det, &                                                          ! determinant
                                     steplength0, & 
                                     steplength, & 
                                     dt, &                                                           ! time increment
                                     aTol
 logical                             error                                                           ! flag indicating an error
 integer(pInt)                       NiterationStress, &                                             ! number of stress integrations
                                     ierr, &                                                         ! error indicator for LAPACK
                                     n, &
                                     o, &
                                     p, &
                                     jacoCounter                                                     ! counter to check for Jacobian update
 integer(pLongInt)                   tick, &
                                     tock, &
                                     tickrate, &
                                     maxticks
  external :: dgesv
  
 !* be pessimistic
 
 crystallite_integrateStress = .false.
#ifndef _OPENMP
 if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
     .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
            .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
   write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> integrateStress at el ip g ',e,i,g
 endif
#endif
 
 
 !* only integrate over fraction of timestep?
 
 if (present(timeFraction)) then
   dt = crystallite_subdt(g,i,e) * timeFraction
   Fg_new = crystallite_subF0(1:3,1:3,g,i,e) + (crystallite_subF(1:3,1:3,g,i,e) - crystallite_subF0(1:3,1:3,g,i,e)) * timeFraction
 else     
   dt = crystallite_subdt(g,i,e)
   Fg_new = crystallite_subF(1:3,1:3,g,i,e)
 endif
 
  
 !* feed local variables
 
 Fp_current =   crystallite_subFp0(1:3,1:3,g,i,e)                   ! "Fp_current" is only used as temp var here...
 Lpguess_old =  crystallite_Lp(1:3,1:3,g,i,e)                       ! consider present Lp good (i.e. worth remembering) ...
 Lpguess =      crystallite_Lp(1:3,1:3,g,i,e)                       ! ... and take it as first guess
 
 
 !* inversion of Fp_current...
 
 invFp_current = math_inv33(Fp_current)                            
 if (all(invFp_current == 0.0_pReal)) then                          ! ... failed?
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,i2,1x,i3)') '<< CRYST >> integrateStress failed on inversion of Fp_current at el ip g ',e,i,g
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) > 0_pInt) then
       write(6,*)
       write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp_current',math_transpose33(Fp_current(1:3,1:3))
     endif
   endif
#endif
   return
 endif
 A = math_mul33x33(Fg_new,invFp_current)                                    ! intermediate tensor needed later to calculate dFe_dLp
  
 
 !* start LpLoop with normal step length
 
 NiterationStress = 0_pInt
 jacoCounter = 0_pInt
 steplength0 = 1.0_pReal
 steplength = steplength0
 residuum_old = 0.0_pReal
 
 LpLoop: do
   NiterationStress = NiterationStress + 1_pInt
    
   
   !* too many loops required ?
   
   if (NiterationStress > nStress) then
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then 
       write(6,'(a,i3,a,i8,1x,i2,1x,i3)') '<< CRYST >> integrateStress reached loop limit',nStress,' at el ip g ',e,i,g
       write(6,*)
     endif
#endif
     return
   endif
    
 
   !* calculate (elastic) 2nd Piola--Kirchhoff stress tensor and its tangent from constitutive law
 
   B = math_I3 - dt*Lpguess
   Fe = math_mul33x33(A,B)                                                    ! current elastic deformation tensor
   call constitutive_TandItsTangent(Tstar, dT_dFe3333, Fe, g,i,e)             ! call constitutive law to calculate 2nd Piola-Kirchhoff stress and its derivative
   Tstar_v = math_Mandel33to6(Tstar)
   p_hydro = sum(Tstar_v(1:3)) / 3.0_pReal
   forall(n=1_pInt:3_pInt) Tstar_v(n) = Tstar_v(n) - p_hydro                  ! get deviatoric stress tensor
    
   
   !* calculate plastic velocity gradient and its tangent from constitutive law
   
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) &
     call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
 
   call constitutive_LpAndItsTangent(Lp_constitutive, dLp_dT_constitutive, Tstar_v, crystallite_Temperature(g,i,e), g, i, e)
 
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
     !$OMP CRITICAL (debugTimingLpTangent)
       debug_cumLpCalls = debug_cumLpCalls + 1_pInt
       debug_cumLpTicks = debug_cumLpTicks + tock-tick
       !$OMP FLUSH (debug_cumLpTicks)
       if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks
     !$OMP END CRITICAL (debugTimingLpTangent)
   endif
    
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
       .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
              .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
     write(6,'(a,i3)') '<< CRYST >> iteration ', NiterationStress
     write(6,*)
     write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lp_constitutive', math_transpose33(Lp_constitutive)
     write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lpguess', math_transpose33(Lpguess)
   endif
#endif
 
 
   !* update current residuum and check for convergence of loop
   
   aTol = max(rTol_crystalliteStress * max(math_norm33(Lpguess),math_norm33(Lp_constitutive)), &    ! absolute tolerance from largest acceptable relative error
              aTol_crystalliteStress)                                                               ! minimum lower cutoff
   residuum = Lpguess - Lp_constitutive
 
   if (any(residuum /= residuum)) then                                                              ! NaN in residuum...
#ifndef _OPENMP
     if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then 
       write(6,'(a,i8,1x,i2,1x,i3,a,i3,a)') '<< CRYST >> integrateStress encountered NaN at el ip g ',e,i,g,&
                                     ' ; iteration ', NiterationStress,&
                                     ' >> returning..!'
     endif
#endif
     return                                                                                         ! ...me = .false. to inform integrator about problem
   elseif (math_norm33(residuum) < aTol) then                                                       ! converged if below absolute tolerance
     exit LpLoop                                                                                    ! ...leave iteration loop
   elseif (math_norm33(residuum) < math_norm33(residuum_old) .or. NiterationStress == 1_pInt ) then ! not converged, but improved norm of residuum (always proceed in first iteration)...
     residuum_old = residuum                                                                        ! ...remember old values and...
     Lpguess_old = Lpguess 
     steplength = steplength0                                                                       ! ...proceed with normal step length (calculate new search direction)
   else                                                                                             ! not converged and residuum not improved...
     steplength = 0.5_pReal * steplength                                                            ! ...try with smaller step length in same direction
     Lpguess = Lpguess_old + steplength * deltaLp
     cycle LpLoop
   endif
 
 
   !* calculate Jacobian for correction term 
   
   if (mod(jacoCounter, iJacoLpresiduum) == 0_pInt) then
     dFe_dLp3333 = 0.0_pReal
     do o=1_pInt,3_pInt; do p=1_pInt,3_pInt
       dFe_dLp3333(p,o,1:3,p) = A(o,1:3)                                                            ! dFe_dLp(i,j,k,l) = -dt * A(i,k) delta(j,l)
     enddo; enddo
     dFe_dLp3333 = -dt * dFe_dLp3333
     dFe_dLp = math_Plain3333to99(dFe_dLp3333)
     dT_dFe_constitutive = math_Plain3333to99(dT_dFe3333)
     dR_dLp = math_identity2nd(9_pInt) - &
              math_mul99x99(dLp_dT_constitutive, math_mul99x99(dT_dFe_constitutive , dFe_dLp))    
     dR_dLp2 = dR_dLp                                                                               ! will be overwritten in first call to LAPACK routine
     work = math_plain33to9(residuum)
#if(FLOAT==8)
     call dgesv(9,1,dR_dLp2,9,ipiv,work,9,ierr)                                                     ! solve dR/dLp * delta Lp = -res for dR/dLp
#elif(FLOAT==4)
     call sgesv(9,1,dR_dLp2,9,ipiv,work,9,ierr)                                                     ! solve dR/dLp * delta Lp = -res for dR/dLp
#endif
     if (ierr /= 0_pInt) then
#ifndef _OPENMP
       if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
         write(6,'(a,i8,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on dR/dLp inversion at el ip g ',e,i,g
         if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
             .and. ((e == debug_e .and. i == debug_i .and. g == debug_g)&
                    .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
           write(6,*)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dR_dLp',transpose(dR_dLp)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dFe_dLp',transpose(dFe_dLp)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dT_dFe_constitutive',transpose(dT_dFe_constitutive)
           write(6,'(a,/,9(12x,9(e15.3,1x)/))') '<< CRYST >> dLp_dT_constitutive',transpose(dLp_dT_constitutive)
           write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> A',math_transpose33(A)
           write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> B',math_transpose33(B)
           write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lp_constitutive',math_transpose33(Lp_constitutive)
           write(6,'(a,/,3(12x,3(e20.7,1x)/))') '<< CRYST >> Lpguess',math_transpose33(Lpguess)
         endif
       endif
#endif
       return
     endif
     deltaLp = - math_plain9to33(work)
   endif
   jacoCounter = jacoCounter + 1_pInt                             ! increase counter for jaco update
 
   Lpguess = Lpguess + steplength * deltaLp
 
 enddo LpLoop
 
 
 !* calculate new plastic and elastic deformation gradient
 
 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new/math_det33(invFp_new)**(1.0_pReal/3.0_pReal)  ! regularize by det
 call math_invert33(invFp_new,Fp_new,det,error)
 if (error .or. any(Fp_new/=Fp_new)) then
#ifndef _OPENMP
   if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
     write(6,'(a,i8,1x,i2,1x,i3,a,i3)') '<< CRYST >> integrateStress failed on invFp_new inversion at el ip g ',&
                                                      e,i,g, ' ; iteration ', NiterationStress
     if (iand(debug_level(debug_crystallite), debug_levelExtensive) /= 0_pInt &
         .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
                .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then
       write(6,*)
       write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> invFp_new',math_transpose33(invFp_new)
     endif
   endif
#endif
   return
 endif
 Fe_new = math_mul33x33(Fg_new,invFp_new)                             ! calc resulting Fe
 
 
 !* add volumetric component to 2nd Piola-Kirchhoff stress and calculate 1st Piola-Kirchhoff stress
 
 forall (n=1_pInt:3_pInt) Tstar_v(n) = Tstar_v(n) + p_hydro
 crystallite_P(1:3,1:3,g,i,e) = math_mul33x33(Fe_new, math_mul33x33(math_Mandel6to33(Tstar_v), math_transpose33(invFp_new)))
  
 
 !* store local values in global variables
 
 crystallite_Lp(1:3,1:3,g,i,e) = Lpguess
 crystallite_Tstar_v(1:6,g,i,e) = Tstar_v
 crystallite_Fp(1:3,1:3,g,i,e) = Fp_new
 crystallite_Fe(1:3,1:3,g,i,e) = Fe_new
 crystallite_invFp(1:3,1:3,g,i,e) = invFp_new
 
 
 !* set return flag to true
 
 crystallite_integrateStress = .true.
#ifndef _OPENMP
 if (iand(debug_level(debug_crystallite),debug_levelExtensive) /= 0_pInt &
     .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) &
             .or. .not. iand(debug_level(debug_crystallite), debug_levelSelective) /= 0_pInt)) then 
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> P / MPa',math_transpose33(crystallite_P(1:3,1:3,g,i,e))/1.0e6_pReal
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Cauchy / MPa', &
              math_mul33x33(crystallite_P(1:3,1:3,g,i,e), math_transpose33(Fg_new)) / 1.0e6_pReal / math_det33(Fg_new)
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fe Lp Fe^-1', &
              math_transpose33(math_mul33x33(Fe_new, math_mul33x33(crystallite_Lp(1:3,1:3,g,i,e), math_inv33(Fe_new))))    ! transpose to get correct print out order
   write(6,'(a,/,3(12x,3(f12.7,1x)/))') '<< CRYST >> Fp',math_transpose33(crystallite_Fp(1:3,1:3,g,i,e))
 endif
#endif
 
 if (iand(debug_level(debug_crystallite), debug_levelBasic) /= 0_pInt) then
   !$OMP CRITICAL (distributionStress)
    debug_StressLoopDistribution(NiterationStress,numerics_integrationMode) = &
      debug_StressLoopDistribution(NiterationStress,numerics_integrationMode) + 1_pInt
   !$OMP END CRITICAL (distributionStress)
 endif
 
end function crystallite_integrateStress
 
 
!--------------------------------------------------------------------------------------------------
!> @brief calculates orientations and disorientations (in case of single grain ips)
!--------------------------------------------------------------------------------------------------
subroutine crystallite_orientations
 use math, only:                       math_pDecomposition, &
                                       math_RtoQ, &
                                       math_qDisorientation, &
                                       math_qConj
 use FEsolving, only:                  FEsolving_execElem, & 
                                       FEsolving_execIP
 use IO, only:                         IO_warning
 use material, only:                   material_phase, &
                                       homogenization_Ngrains, &
                                       phase_localPlasticity, &
                                       phase_plasticityInstance
 use mesh, only:                       mesh_element, &
                                       mesh_ipNeighborhood, &
                                       FE_NipNeighbors, &
                                       FE_geomtype
 use constitutive_nonlocal, only:      constitutive_nonlocal_structure, &
                                       constitutive_nonlocal_updateCompatibility
 
 implicit none
 integer(pInt)                   e, &                          ! element index
                                 i, &                          ! integration point index
                                 g, &                          ! grain index
                                 n, &                          ! neighbor index 
                                 neighboring_e, &              ! element index of my neighbor
                                 neighboring_i, &              ! integration point index of my neighbor
                                 myPhase, &                    ! phase
                                 neighboringPhase, &
                                 myInstance, &                 ! instance of plasticity
                                 neighboringInstance, &
                                 myStructure, &                ! lattice structure
                                 neighboringStructure
 real(pReal), dimension(3,3) ::  U, R
 real(pReal), dimension(4) ::    orientation
 logical error
 
 ! --- CALCULATE ORIENTATION AND LATTICE ROTATION ---
 
 !$OMP PARALLEL DO PRIVATE(error,U,R,orientation)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       do g = 1_pInt,homogenization_Ngrains(mesh_element(3,e))
         
         call math_pDecomposition(crystallite_Fe(1:3,1:3,g,i,e), U, R, error)                              ! polar decomposition of Fe
         if (error) then
           call IO_warning(650_pInt, e, i, g)
           orientation = [1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]                                ! fake orientation
         else
           orientation = math_RtoQ(transpose(R))
         endif
         crystallite_rotation(1:4,g,i,e) = math_qDisorientation(crystallite_orientation0(1:4,g,i,e), &  ! active rotation from ori0
                                                                         orientation, &                          ! to current orientation
                                                                         0_pInt )                                ! we don't want symmetry here  
         crystallite_orientation(1:4,g,i,e) = orientation
       enddo
     enddo
   enddo
 !$OMP END PARALLEL DO
 
 
 ! --- UPDATE SOME ADDITIONAL VARIABLES THAT ARE NEEDED FOR NONLOCAL MATERIAL ---
 ! --- we use crystallite_orientation from above, so need a seperate loop
 
 !$OMP PARALLEL DO PRIVATE(myPhase,myInstance,myStructure,neighboring_e,neighboring_i,neighboringPhase,&
 !$OMP neighboringInstance,neighboringStructure)
   do e = FEsolving_execElem(1),FEsolving_execElem(2)
     do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
       myPhase = material_phase(1,i,e)                                                                     ! get my phase
       if (.not. phase_localPlasticity(myPhase)) then                                                      ! if nonlocal model
         myInstance = phase_plasticityInstance(myPhase)
         myStructure = constitutive_nonlocal_structure(myInstance)                                         ! get my crystal structure
 
         
         ! --- calculate disorientation between me and my neighbor ---
         
         do n = 1_pInt,FE_NipNeighbors(FE_geomtype(mesh_element(2,e)))                                     ! loop through my neighbors
           neighboring_e = mesh_ipNeighborhood(1,n,i,e)
           neighboring_i = mesh_ipNeighborhood(2,n,i,e)
           if ((neighboring_e > 0) .and. (neighboring_i > 0)) then                                         ! if neighbor exists
             neighboringPhase = material_phase(1,neighboring_i,neighboring_e)                              ! get my neighbor's phase
             if (.not. phase_localPlasticity(neighboringPhase)) then                                       ! neighbor got also nonlocal plasticity
               neighboringInstance = phase_plasticityInstance(neighboringPhase)        
               neighboringStructure = constitutive_nonlocal_structure(neighboringInstance)                 ! get my neighbor's crystal structure
               if (myStructure == neighboringStructure) then                                               ! if my neighbor has same crystal structure like me
                 crystallite_disorientation(:,n,1,i,e) = &
                   math_qDisorientation( crystallite_orientation(1:4,1,i,e), &
                                                  crystallite_orientation(1:4,1,neighboring_i,neighboring_e), & 
                                                  crystallite_symmetryID(1,i,e))                           ! calculate disorientation
               else                                                                                        ! for neighbor with different phase
                 crystallite_disorientation(:,n,1,i,e) = (/0.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal/)    ! 180 degree rotation about 100 axis
               endif
             else                                                                                          ! for neighbor with local plasticity
               crystallite_disorientation(:,n,1,i,e) = (/-1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal/)     ! homomorphic identity
             endif
           else                                                                                            ! no existing neighbor
             crystallite_disorientation(:,n,1,i,e) = (/-1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal/)       ! homomorphic identity
           endif
         enddo
 
 
         ! --- calculate compatibility and transmissivity between me and my neighbor ---
 
         call constitutive_nonlocal_updateCompatibility(crystallite_orientation,i,e)
 
       endif
     enddo
   enddo
 !$OMP END PARALLEL DO
 
end subroutine crystallite_orientations


 
!--------------------------------------------------------------------------------------------------
!> @brief return results of particular grain
!--------------------------------------------------------------------------------------------------
function crystallite_postResults(&
   dt,&             ! time increment
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )
 use math, only:                      math_qToEuler, &
                                      math_qToAxisAngle, &
                                      math_mul33x33, &
                                      math_transpose33, &
                                      math_det33, &
                                      math_I3, &
                                      inDeg, &
                                      math_Mandel6to33, &
                                      math_qMul, &
                                      math_qConj
 use mesh, only:                      mesh_element, &
                                      mesh_ipVolume
 use material, only:                  microstructure_crystallite, &
                                      crystallite_Noutput, &
                                      material_phase, &
                                      material_texture, &
                                      homogenization_Ngrains
 use constitutive, only:              constitutive_sizePostResults, &
                                      constitutive_postResults, &
                                      constitutive_homogenizedC
 
 implicit none
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 real(pReal), intent(in)::            dt                            ! time increment

 !*** output variables ***!
 real(pReal), dimension(1+crystallite_sizePostResults(microstructure_crystallite(mesh_element(4,e)))+ &
                        1+constitutive_sizePostResults(g,i,e)) :: crystallite_postResults
 
 !*** local variables ***!
 real(pReal), dimension(3,3) ::       Ee
 real(pReal), dimension(4) ::         rotation
 real(pReal)                          detF
 integer(pInt)                        o,c,crystID,mySize

 crystID = microstructure_crystallite(mesh_element(4,e))

 crystallite_postResults = 0.0_pReal
 c = 0_pInt
 crystallite_postResults(c+1) = real(crystallite_sizePostResults(crystID),pReal)                                  ! size of results from cryst
 c = c + 1_pInt
 
 do o = 1_pInt,crystallite_Noutput(crystID)
   mySize = 0_pInt
   select case(crystallite_output(o,crystID))
     case ('phase')
       mySize = 1_pInt
       crystallite_postResults(c+1) = real(material_phase(g,i,e),pReal)                                           ! phaseID of grain
     case ('texture')
       mySize = 1_pInt
       crystallite_postResults(c+1) = real(material_texture(g,i,e),pReal)                                         ! textureID of grain
     case ('volume')
       mySize = 1_pInt
       detF = math_det33(crystallite_partionedF(1:3,1:3,g,i,e))                                                   ! V_current = det(F) * V_reference
       crystallite_postResults(c+1) = detF * mesh_ipVolume(i,e) / homogenization_Ngrains(mesh_element(3,e))       ! grain volume (not fraction but absolute)
     case ('orientation')
       mySize = 4_pInt
       crystallite_postResults(c+1:c+mySize) = crystallite_orientation(1:4,g,i,e)                                 ! grain orientation as quaternion
     case ('eulerangles')
       mySize = 3_pInt
       crystallite_postResults(c+1:c+mySize) = inDeg * math_qToEuler(crystallite_orientation(1:4,g,i,e)) ! grain orientation as Euler angles in degree
     case ('grainrotation')
       mySize = 4_pInt
       crystallite_postResults(c+1:c+mySize) = math_qToAxisAngle(crystallite_rotation(1:4,g,i,e))        ! grain rotation away from initial orientation as axis-angle in crystal reference coordinates
       crystallite_postResults(c+4) = inDeg * crystallite_postResults(c+4)                                        ! angle in degree
     case ('grainrotationx')
       mySize = 1_pInt
       rotation = math_qToAxisAngle(math_qMul(math_qMul(crystallite_orientation(1:4,g,i,e), &
                                                                 crystallite_rotation(1:4,g,i,e)), &
                                                                 math_qConj(crystallite_orientation(1:4,g,i,e)))) ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+1) = inDeg * rotation(1) * rotation(4)                                           ! angle in degree
     case ('grainrotationy')
       mySize = 1_pInt
       rotation = math_qToAxisAngle(math_qMul(math_qMul(crystallite_orientation(1:4,g,i,e), &
                                                                 crystallite_rotation(1:4,g,i,e)), &
                                                                 math_qConj(crystallite_orientation(1:4,g,i,e)))) ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+1) = inDeg * rotation(2) * rotation(4)                                           ! angle in degree
     case ('grainrotationz')
       mySize = 1_pInt
       rotation = math_qToAxisAngle(math_qMul(math_qMul(crystallite_orientation(1:4,g,i,e), &
                                                                 crystallite_rotation(1:4,g,i,e)), &
                                                                 math_qConj(crystallite_orientation(1:4,g,i,e)))) ! grain rotation away from initial orientation as axis-angle in sample reference coordinates
       crystallite_postResults(c+1) = inDeg * rotation(3) * rotation(4)                                           ! angle in degree

! remark: tensor output is of the form 11,12,13, 21,22,23, 31,32,33
! thus row index i is slow, while column index j is fast. reminder: "row is slow"
  
     case ('defgrad','f')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(math_transpose33(crystallite_partionedF(1:3,1:3,g,i,e)),[mySize])
     case ('e')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = 0.5_pReal * reshape((math_mul33x33( &
                                               math_transpose33(crystallite_partionedF(1:3,1:3,g,i,e)), &
                                               crystallite_partionedF(1:3,1:3,g,i,e)) - math_I3),[mySize])
     case ('fe')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(math_transpose33(crystallite_Fe(1:3,1:3,g,i,e)),[mySize])
     case ('ee')
       Ee = 0.5_pReal * (math_mul33x33(math_transpose33(crystallite_Fe(1:3,1:3,g,i,e)), crystallite_Fe(1:3,1:3,g,i,e)) - math_I3)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(Ee,[mySize])
     case ('fp')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(math_transpose33(crystallite_Fp(1:3,1:3,g,i,e)),[mySize])
     case ('lp')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(math_transpose33(crystallite_Lp(1:3,1:3,g,i,e)),[mySize])
     case ('p','firstpiola','1stpiola')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(math_transpose33(crystallite_P(1:3,1:3,g,i,e)),[mySize])
     case ('s','tstar','secondpiola','2ndpiola')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(math_Mandel6to33(crystallite_Tstar_v(1:6,g,i,e)),[mySize])
     case ('elasmatrix')
       mySize = 36_pInt
       crystallite_postResults(c+1:c+mySize) = reshape(constitutive_homogenizedC(g,i,e),(/mySize/))
   end select
   c = c + mySize
 enddo

 crystallite_postResults(c+1) = real(constitutive_sizePostResults(g,i,e),pReal)             ! size of constitutive results
 c = c + 1_pInt
 if (constitutive_sizePostResults(g,i,e) > 0_pInt) &
   crystallite_postResults(c+1:c+constitutive_sizePostResults(g,i,e)) = constitutive_postResults(crystallite_Tstar_v(1:6,g,i,e), &
                                                                                                 crystallite_Fe, &
                                                                                                 crystallite_Temperature(g,i,e), &
                                                                                                 dt, g, i, e)
 c = c + constitutive_sizePostResults(g,i,e)

end function crystallite_postResults

end module crystallite
