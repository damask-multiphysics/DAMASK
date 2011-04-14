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
!***************************************
!*      Module: CRYSTALLITE            *
!***************************************
!* contains:                           *
!* - _init                             *
!* - materialpoint_stressAndItsTangent *
!* - _partitionDeformation             *
!* - _updateState                      *
!* - _stressAndItsTangent              *
!* - _postResults                      *
!***************************************

MODULE crystallite

use prec, only: pReal, pInt
implicit none

! ****************************************************************
! *** General variables for the crystallite calculation        ***
! ****************************************************************

integer(pInt) crystallite_maxSizePostResults
integer(pInt), dimension(:),       allocatable :: crystallite_sizePostResults
integer(pInt), dimension(:,:),     allocatable :: crystallite_sizePostResult
character(len=64), dimension(:,:),   allocatable :: crystallite_output             ! name of each post result output
integer(pInt), dimension (:,:,:),  allocatable :: &    
    crystallite_symmetryID               ! crystallographic symmetry 1=cubic 2=hexagonal, needed in all orientation calcs     
    
real(pReal), dimension (:,:,:), allocatable :: &
    crystallite_dt, &                    ! requested time increment of each grain
    crystallite_subdt, &                 ! substepped time increment of each grain
    crystallite_subFrac, &               ! already calculated fraction of increment
    crystallite_subStep, &               ! size of next integration step
    crystallite_statedamper, &           ! damping for state update
    crystallite_Temperature, &           ! Temp of each grain
    crystallite_partionedTemperature0, & ! Temp of each grain at start of homog inc
    crystallite_subTemperature0, &       ! Temp of each grain at start of crystallite inc
    crystallite_dotTemperature           ! evolution of Temperature of each grain
real(pReal), dimension (:,:,:,:), allocatable :: &
    crystallite_Tstar_v, &               ! current 2nd Piola-Kirchhoff stress vector (end of converged time step)
    crystallite_Tstar0_v, &              ! 2nd Piola-Kirchhoff stress vector at start of FE inc
    crystallite_partionedTstar0_v, &     ! 2nd Piola-Kirchhoff stress vector at start of homog inc
    crystallite_subTstar0_v, &           ! 2nd Piola-Kirchhoff stress vector at start of crystallite inc
    crystallite_orientation, &           ! orientation as quaternion
    crystallite_orientation0, &          ! initial orientation as quaternion
    crystallite_rotation                 ! grain rotation away from initial orientation as axis-angle (in degrees)
real(pReal), dimension (:,:,:,:,:), allocatable :: &
    crystallite_Fe, &                    ! current "elastic" def grad (end of converged time step)
    crystallite_Fp, &                    ! current plastic def grad (end of converged time step)
    crystallite_invFp, &                 ! inverse of current plastic def grad (end of converged time step)
    crystallite_Fp0, &                   ! plastic def grad at start of FE inc
    crystallite_partionedFp0,&           ! plastic def grad at start of homog inc
    crystallite_subFp0,&                 ! plastic def grad at start of crystallite inc
    crystallite_F0, &                    ! def grad at start of FE inc
    crystallite_partionedF,  &           ! def grad to be reached at end of homog inc
    crystallite_partionedF0, &           ! def grad at start of homog inc
    crystallite_subF,  &                 ! def grad to be reached at end of crystallite inc
    crystallite_subF0, &                 ! def grad at start of crystallite inc
    crystallite_Lp, &                    ! current plastic velocitiy grad (end of converged time step)
    crystallite_Lp0, &                   ! plastic velocitiy grad at start of FE inc
    crystallite_partionedLp0,&           ! plastic velocity grad at start of homog inc
    crystallite_subLp0,&                 ! plastic velocity grad at start of crystallite inc
    crystallite_P, &                     ! 1st Piola-Kirchhoff stress per grain
    crystallite_disorientation           ! disorientation between two neighboring ips (only calculated for single grain IPs)
real(pReal), dimension (:,:,:,:,:,:,:), allocatable :: &
    crystallite_dPdF, &                  ! current individual dPdF per grain (end of converged time step)
    crystallite_dPdF0, &                 ! individual dPdF per grain at start of FE inc
    crystallite_partioneddPdF0, &        ! individual dPdF per grain at start of homog inc
    crystallite_fallbackdPdF             ! dPdF fallback for non-converged grains (elastic prediction)
logical, dimension (:,:,:), allocatable :: &
    crystallite_localConstitution, &     ! indicates this grain to have purely local constitutive law
    crystallite_requested, &             ! flag to request crystallite calculation
    crystallite_todo, &                  ! flag to indicate need for further computation
    crystallite_converged                ! convergence flag

CONTAINS


!********************************************************************
! allocate and initialize per grain variables
!********************************************************************
subroutine crystallite_init(Temperature)
  
!*** variables and functions from other modules ***!
use prec, only:       pInt, &
                      pReal
use debug, only:      debug_info, &
                      debug_reset, &
                      debug_verbosity
use numerics, only:   subStepSizeCryst, &
                      stepIncreaseCryst
use math, only:       math_I3, &
                      math_EulerToR, &
                      math_inv3x3, &
                      math_transpose3x3, &
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
use lattice, only:    lattice_symmetryType, &
                      lattice_Sslip,lattice_Sslip_v,lattice_Stwin,lattice_Stwin_v, lattice_maxNslipFamily, lattice_maxNtwinFamily, &
                      lattice_NslipSystem,lattice_NtwinSystem

use constitutive, only: constitutive_microstructure
use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_label, &
                                      constitutive_phenopowerlaw_structure, &
                                      constitutive_phenopowerlaw_Nslip
use constitutive_titanmod, only:      constitutive_titanmod_label, &
                                      constitutive_titanmod_structure
use constitutive_dislotwin, only:     constitutive_dislotwin_label, &
                                      constitutive_dislotwin_structure
use constitutive_nonlocal, only:      constitutive_nonlocal_label, &
                                      constitutive_nonlocal_structure
 
implicit none
integer(pInt), parameter :: file = 200, &
                            maxNchunks = 2
 
!*** input variables ***!
real(pReal) Temperature

!*** output variables ***!

!*** local variables ***!
integer(pInt), dimension(1+2*maxNchunks) :: positions
integer(pInt)               g, &                          ! grain number
                            i, &                          ! integration point number
                            e, &                          ! element number
                            gMax, &                       ! maximum number of grains
                            iMax, &                       ! maximum number of integration points
                            eMax, &                       ! maximum number of elements
                            nMax, &                       ! maximum number of ip neighbors
                            myNgrains, &                  ! number of grains in current IP
                            section, &
                            j, &
                            p, &
                            output, &
                            mySize, &
                            myStructure, &                ! lattice structure 
                            myPhase, &
                            myMat
character(len=64)           tag
character(len=1024)         line


!$OMP CRITICAL (write2out)
  write(6,*)
  write(6,*) '<<<+-  crystallite init  -+>>>'
  write(6,*) '$Id$'
  write(6,*)
!$OMP END CRITICAL (write2out)


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
allocate(crystallite_statedamper(gMax,iMax,eMax));                     crystallite_statedamper = 1.0_pReal
allocate(crystallite_symmetryID(gMax,iMax,eMax));                       crystallite_symmetryID = 0.0_pReal !NEW
allocate(crystallite_orientation(4,gMax,iMax,eMax));                   crystallite_orientation = 0.0_pReal
allocate(crystallite_orientation0(4,gMax,iMax,eMax));                 crystallite_orientation0 = 0.0_pReal
allocate(crystallite_rotation(4,gMax,iMax,eMax));                         crystallite_rotation = 0.0_pReal
allocate(crystallite_disorientation(4,nMax,gMax,iMax,eMax));        crystallite_disorientation = 0.0_pReal
allocate(crystallite_localConstitution(gMax,iMax,eMax));         crystallite_localConstitution = .true.
allocate(crystallite_requested(gMax,iMax,eMax));                         crystallite_requested = .false.
allocate(crystallite_todo(gMax,iMax,eMax));                                   crystallite_todo = .false.
allocate(crystallite_converged(gMax,iMax,eMax));                         crystallite_converged = .true.
allocate(crystallite_output(maxval(crystallite_Noutput), &
                            material_Ncrystallite)) ;                       crystallite_output = ''
allocate(crystallite_sizePostResults(material_Ncrystallite)) ;     crystallite_sizePostResults = 0_pInt
allocate(crystallite_sizePostResult(maxval(crystallite_Noutput), &
                                    material_Ncrystallite)) ;       crystallite_sizePostResult = 0_pInt


if(.not. IO_open_file(file,material_configFile)) call IO_error (100) ! corrupt config file
line = ''
section = 0

do while (IO_lc(IO_getTag(line,'<','>')) /= material_partCrystallite)     ! wind forward to <crystallite>
  read(file,'(a1024)',END=100) line
enddo

do                                                       ! read thru sections of phase part
  read(file,'(a1024)',END=100) line
  if (IO_isBlank(line)) cycle                            ! skip empty lines
  if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
  if (IO_getTag(line,'[',']') /= '') then                ! next section
    section = section + 1
    output = 0                                           ! reset output counter
  endif
  if (section > 0) then
    positions = IO_stringPos(line,maxNchunks)
    tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
    select case(tag)
      case ('(output)')
        output = output + 1
        crystallite_output(output,section) = IO_lc(IO_stringValue(line,positions,2))
    end select
  endif
enddo

100 close(file)

do i = 1,material_Ncrystallite                        ! sanity checks
enddo

do i = 1,material_Ncrystallite
  do j = 1,crystallite_Noutput(i)
    select case(crystallite_output(j,i))
      case('phase','texture','volume')
        mySize = 1
      case('orientation','grainrotation')   ! orientation as quaternion, or deviation from initial grain orientation in axis-angle form (angle in degrees)
        mySize = 4
      case('eulerangles')   ! Bunge (3-1-3) Euler angles
        mySize = 3
      case('defgrad','f','fe','fp','lp','ee','p','firstpiola','1stpiola','s','tstar','secondpiola','2ndpiola')
        mySize = 9
      case default
        mySize = 0      
    end select
  
    if (mySize > 0_pInt) then                               ! any meaningful output found
      crystallite_sizePostResult(j,i) = mySize
      crystallite_sizePostResults(i) = crystallite_sizePostResults(i) + mySize
    endif
  enddo
enddo
crystallite_maxSizePostResults = maxval(crystallite_sizePostResults)

! write description file for crystallite output

if(.not. IO_open_jobFile(file,'outputCrystallite')) call IO_error (50) ! problems in writing file
 
do p = 1,material_Ncrystallite
  write(file,*)
  write(file,'(a)') '['//trim(crystallite_name(p))//']'
  write(file,*)
  do e = 1,crystallite_Noutput(p)
    write(file,'(a,i4)') trim(crystallite_output(e,p))//char(9),crystallite_sizePostResult(e,p)
  enddo
enddo

close(file)

!$OMP PARALLEL PRIVATE(myNgrains,myPhase,myMat,myStructure)

!$OMP DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                      ! iterate over all cp elements
    myNgrains = homogenization_Ngrains(mesh_element(3,e))                                                 ! look up homogenization-->grainCount
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                    ! iterate over IPs of this element
      do g = 1,myNgrains
        crystallite_Fp0(1:3,1:3,g,i,e) = math_EulerToR(material_EulerAngles(1:3,g,i,e))                   ! plastic def gradient reflects init orientation
        crystallite_F0(1:3,1:3,g,i,e)  = math_I3
        crystallite_localConstitution(g,i,e) = phase_localConstitution(material_phase(g,i,e))
        !$OMP FLUSH(crystallite_Fp0)
        crystallite_Fe(1:3,1:3,g,i,e)  = math_transpose3x3(crystallite_Fp0(1:3,1:3,g,i,e))
      enddo
    enddo
  enddo
!$OMP ENDDO
crystallite_partionedTemperature0 = Temperature   ! isothermal assumption
crystallite_partionedFp0 = crystallite_Fp0
crystallite_partionedF0 = crystallite_F0
crystallite_partionedF = crystallite_F0
crystallite_requested = .true.


! Initialize crystallite_symmetryID(g,i,e)

!$OMP DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do g = 1,myNgrains
        myPhase = material_phase(g,i,e)
        myMat   = phase_constitutionInstance(myPhase)
        select case (phase_constitution(myPhase))
          case (constitutive_phenopowerlaw_label)
            myStructure = constitutive_phenopowerlaw_structure(myMat)
          case (constitutive_titanmod_label)
            myStructure = constitutive_titanmod_structure(myMat)
          case (constitutive_dislotwin_label)
            myStructure = constitutive_dislotwin_structure(myMat)
          case (constitutive_nonlocal_label)
            myStructure = constitutive_nonlocal_structure(myMat)
          case default
            myStructure = -1_pInt ! does this happen for j2 material?
        end select
        if (myStructure > 0_pInt) then   
          crystallite_symmetryID(g,i,e) = lattice_symmetryType(myStructure) ! structure = 1(fcc) or 2(bcc) => 1; 3(hex)=>2  
        endif
      enddo
    enddo
  enddo
!$OMP ENDDO   

!$OMP END PARALLEL
 
call crystallite_orientations()
crystallite_orientation0 = crystallite_orientation             ! Store initial orientations for calculation of grain rotations

!$OMP PARALLEL DO PRIVATE(myNgrains)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do g = 1,myNgrains
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)  ! update dependent state variables to be consistent with basic states    
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO

call crystallite_stressAndItsTangent(.true.)                   ! request elastic answers
crystallite_fallbackdPdF = crystallite_dPdF                    ! use initial elastic stiffness as fallback
 
!    *** Output to MARC output file ***
if (debug_verbosity > 0) then
  !$OMP CRITICAL (write2out)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Temperature:           ', shape(crystallite_Temperature)
    write(6,'(a35,x,7(i5,x))') 'crystallite_dotTemperature:        ', shape(crystallite_dotTemperature)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Fe:                    ', shape(crystallite_Fe)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Fp:                    ', shape(crystallite_Fp)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Lp:                    ', shape(crystallite_Lp)
    write(6,'(a35,x,7(i5,x))') 'crystallite_F0:                    ', shape(crystallite_F0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Fp0:                   ', shape(crystallite_Fp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Lp0:                   ', shape(crystallite_Lp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedF:            ', shape(crystallite_partionedF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedTemp0:        ', shape(crystallite_partionedTemperature0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedF0:           ', shape(crystallite_partionedF0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedFp0:          ', shape(crystallite_partionedFp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedLp0:          ', shape(crystallite_partionedLp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subF:                  ', shape(crystallite_subF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subTemperature0:       ', shape(crystallite_subTemperature0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_symmetryID:            ', shape(crystallite_symmetryID)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subF0:                 ', shape(crystallite_subF0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subFp0:                ', shape(crystallite_subFp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subLp0:                ', shape(crystallite_subLp0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_P:                     ', shape(crystallite_P)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Tstar_v:               ', shape(crystallite_Tstar_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_Tstar0_v:              ', shape(crystallite_Tstar0_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partionedTstar0_v:     ', shape(crystallite_partionedTstar0_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subTstar0_v:           ', shape(crystallite_subTstar0_v)
    write(6,'(a35,x,7(i5,x))') 'crystallite_dPdF:                  ', shape(crystallite_dPdF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_dPdF0:                 ', shape(crystallite_dPdF0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_partioneddPdF0:        ', shape(crystallite_partioneddPdF0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_fallbackdPdF:          ', shape(crystallite_fallbackdPdF)
    write(6,'(a35,x,7(i5,x))') 'crystallite_orientation:           ', shape(crystallite_orientation)
    write(6,'(a35,x,7(i5,x))') 'crystallite_orientation0:          ', shape(crystallite_orientation0)
    write(6,'(a35,x,7(i5,x))') 'crystallite_rotation:              ', shape(crystallite_rotation)
    write(6,'(a35,x,7(i5,x))') 'crystallite_disorientation:        ', shape(crystallite_disorientation)
    write(6,'(a35,x,7(i5,x))') 'crystallite_dt:                    ', shape(crystallite_dt)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subdt:                 ', shape(crystallite_subdt)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subFrac:               ', shape(crystallite_subFrac)
    write(6,'(a35,x,7(i5,x))') 'crystallite_subStep:               ', shape(crystallite_subStep)
    write(6,'(a35,x,7(i5,x))') 'crystallite_stateDamper:           ', shape(crystallite_stateDamper)
    write(6,'(a35,x,7(i5,x))') 'crystallite_localConstitution:     ', shape(crystallite_localConstitution)
    write(6,'(a35,x,7(i5,x))') 'crystallite_requested:             ', shape(crystallite_requested)
    write(6,'(a35,x,7(i5,x))') 'crystallite_todo:                  ', shape(crystallite_todo)
    write(6,'(a35,x,7(i5,x))') 'crystallite_converged:             ', shape(crystallite_converged)
    write(6,'(a35,x,7(i5,x))') 'crystallite_sizePostResults:       ', shape(crystallite_sizePostResults)
    write(6,'(a35,x,7(i5,x))') 'crystallite_sizePostResult:        ', shape(crystallite_sizePostResult)
    write(6,*)
    write(6,*) 'Number of nonlocal grains: ',count(.not. crystallite_localConstitution)
    call flush(6)
  !$OMP END CRITICAL (write2out)
endif

call debug_info()
call debug_reset()

endsubroutine


 
!********************************************************************
! calculate stress (P) and tangent (dPdF) for crystallites
!********************************************************************
subroutine crystallite_stressAndItsTangent(updateJaco)

!*** variables and functions from other modules ***!
use prec, only:                                       pInt, &
                                                      pReal
use numerics, only:                                   subStepMinCryst, &
                                                      subStepSizeCryst, &
                                                      stepIncreaseCryst, &
                                                      pert_Fg, &
                                                      pert_method, &
                                                      nCryst, &
                                                      iJacoStiffness, &
                                                      numerics_integrator, &
                                                      numerics_integrationMode
use debug, only:                                      debug_verbosity, &
                                                      debug_selectiveDebugger, &
                                                      debug_e, &
                                                      debug_i, &
                                                      debug_g, &
                                                      debug_CrystalliteLoopDistribution
use IO, only:                                         IO_warning
use math, only:                                       math_inv3x3, &
                                                      math_transpose3x3, &
                                                      math_mul33x33, &
                                                      math_mul66x6, &
                                                      math_Mandel6to33, &
                                                      math_Mandel33to6, &
                                                      math_transpose3x3, &
                                                      math_I3
use FEsolving, only:                                  FEsolving_execElem, & 
                                                      FEsolving_execIP
use mesh, only:                                       mesh_element, &
                                                      mesh_NcpElems, &
                                                      mesh_maxNips
use material, only:                                   homogenization_Ngrains, &
                                                      homogenization_maxNgrains
use constitutive, only:                               constitutive_sizeState, &
                                                      constitutive_sizeDotState, &
                                                      constitutive_state, &
                                                      constitutive_state_backup, &
                                                      constitutive_subState0, &
                                                      constitutive_partionedState0, &
                                                      constitutive_homogenizedC, &
                                                      constitutive_dotState, &
                                                      constitutive_dotState_backup

implicit none

!*** input variables ***!
logical, intent(in) ::                                updateJaco                    ! flag indicating wehther we want to update the Jacobian (stiffness) or not

!*** output variables ***!

!*** local variables ***!
real(pReal)                                           myPert, &                     ! perturbation with correct sign
                                                      formerSubStep
real(pReal), dimension(3,3) ::                        invFp, &                      ! inverse of the plastic deformation gradient
                                                      Fe_guess, &                   ! guess for elastic deformation gradient
                                                      Tstar                         ! 2nd Piola-Kirchhoff stress tensor
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
integer(pInt)                                         NiterationCrystallite, &      ! number of iterations in crystallite loop
                                                      e, &                          ! element index
                                                      i, &                          ! integration point index
                                                      g, &                          ! grain index
                                                      k, &
                                                      l, &
                                                      perturbation , &              ! loop counter for forward,backward perturbation mode
                                                      myNgrains, &
                                                      mySizeState, &
                                                      mySizeDotState
logical, dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                                      convergenceFlag_backup

! --+>> INITIALIZE TO STARTING CONDITION <<+--

if (debug_verbosity > 4 .and. debug_e > 0 .and. debug_e <= mesh_NcpElems &
                        .and. debug_i > 0 .and. debug_i <= mesh_maxNips &
                        .and. debug_g > 0 .and. debug_g <= homogenization_maxNgrains) then
  !$OMP CRITICAL (write2out)
    write (6,*)
    write (6,'(a,i5,x,i2,x,i3)') '<< CRYST >> crystallite start at el ip g ', debug_e, debug_i, debug_g
    write (6,'(a,/,12(x),f14.9)') '<< CRYST >> Temp0', crystallite_partionedTemperature0(debug_g,debug_i,debug_e)
    write (6,'(a,/,3(12(x),3(f14.9,x)/))') '<< CRYST >> F0 ', &
                                          math_transpose3x3(crystallite_partionedF0(1:3,1:3,debug_g,debug_i,debug_e))
    write (6,'(a,/,3(12(x),3(f14.9,x)/))') '<< CRYST >> Fp0', &
                                          math_transpose3x3(crystallite_partionedFp0(1:3,1:3,debug_g,debug_i,debug_e))
    write (6,'(a,/,3(12(x),3(f14.9,x)/))') '<< CRYST >> Lp0', &
                                          math_transpose3x3(crystallite_partionedLp0(1:3,1:3,debug_g,debug_i,debug_e))
  !$OMP END CRITICAL (write2out)
endif

crystallite_subStep = 0.0_pReal

!$OMP PARALLEL DO PRIVATE(myNgrains)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
      do g = 1,myNgrains
        if (crystallite_requested(g,i,e)) then                                                          ! initialize restoration point of ...
          crystallite_subTemperature0(g,i,e) = crystallite_partionedTemperature0(g,i,e)                 ! ...temperature
          constitutive_subState0(g,i,e)%p = constitutive_partionedState0(g,i,e)%p                       ! ...microstructure
          crystallite_subFp0(1:3,1:3,g,i,e) = crystallite_partionedFp0(1:3,1:3,g,i,e)                   ! ...plastic def grad
          crystallite_subLp0(1:3,1:3,g,i,e) = crystallite_partionedLp0(1:3,1:3,g,i,e)                   ! ...plastic velocity grad
          crystallite_dPdF0(1:3,1:3,1:3,1:3,g,i,e) = crystallite_partioneddPdF0(1:3,1:3,1:3,1:3,g,i,e)  ! ...stiffness
          crystallite_subF0(1:3,1:3,g,i,e) = crystallite_partionedF0(1:3,1:3,g,i,e)                     ! ...def grad
          crystallite_subTstar0_v(1:6,g,i,e) = crystallite_partionedTstar0_v(1:6,g,i,e)                 !...2nd PK stress

          crystallite_subFrac(g,i,e) = 0.0_pReal
          crystallite_subStep(g,i,e) = 1.0_pReal/subStepSizeCryst
          crystallite_todo(g,i,e) = .true.
          crystallite_converged(g,i,e) = .false.                                                        ! pretend failed step of twice the required size
        endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO


! --+>> CRYSTALLITE CUTBACK LOOP <<+--

NiterationCrystallite = 0_pInt
numerics_integrationMode = 1_pInt
do while (any(crystallite_subStep(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinCryst))  ! cutback loop for crystallites

  !$OMP PARALLEL DO PRIVATE(myNgrains,formerSubStep)
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                  ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          
          ! --- wind forward ---
          
          if (crystallite_converged(g,i,e)) then
#ifndef _OPENMP
            if (debug_verbosity > 4 &
                .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
              write(6,'(a,f10.8,a,f10.8,a)') '<< CRYST >> winding forward from ', &
                crystallite_subFrac(g,i,e),' to current crystallite_subfrac ', &
                crystallite_subFrac(g,i,e)+crystallite_subStep(g,i,e),' in crystallite_stressAndItsTangent'
              write(6,*)
            endif
#endif
            crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
            formerSubStep = crystallite_subStep(g,i,e)
            !$OMP FLUSH(crystallite_subFrac)
            crystallite_subStep(g,i,e) = min( 1.0_pReal - crystallite_subFrac(g,i,e), &
                                              stepIncreaseCryst * crystallite_subStep(g,i,e) )
            !$OMP FLUSH(crystallite_subStep)
            if (crystallite_subStep(g,i,e) > subStepMinCryst) then
              crystallite_subTemperature0(g,i,e) = crystallite_Temperature(g,i,e)                       ! wind forward...
              crystallite_subF0(1:3,1:3,g,i,e) = crystallite_subF(1:3,1:3,g,i,e)                        ! ...def grad
              crystallite_subFp0(1:3,1:3,g,i,e) = crystallite_Fp(1:3,1:3,g,i,e)                         ! ...plastic def grad
              crystallite_subLp0(1:3,1:3,g,i,e) = crystallite_Lp(1:3,1:3,g,i,e)                         ! ...plastic velocity gradient
              constitutive_subState0(g,i,e)%p = constitutive_state(g,i,e)%p                             ! ...microstructure
              crystallite_subTstar0_v(1:6,g,i,e) = crystallite_Tstar_v(1:6,g,i,e)                       ! ...2nd PK stress
              !$OMP FLUSH(crystallite_subF0)
            elseif (formerSubStep > subStepMinCryst) then                                               ! this crystallite just converged
              if (debug_verbosity > 0) then
                !$OMP CRITICAL (distributionCrystallite)
                  debug_CrystalliteLoopDistribution(min(nCryst+1,NiterationCrystallite)) = &
                    debug_CrystalliteLoopDistribution(min(nCryst+1,NiterationCrystallite)) + 1
                !$OMP END CRITICAL (distributionCrystallite)
              endif
            endif
          
          ! --- cutback ---
          
          else
            crystallite_subStep(g,i,e) = subStepSizeCryst * crystallite_subStep(g,i,e)                  ! cut step in half and restore...
            crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e)                         ! ...temperature
            crystallite_Fp(1:3,1:3,g,i,e) = crystallite_subFp0(1:3,1:3,g,i,e)                           ! ...plastic def grad
            crystallite_invFp(1:3,1:3,g,i,e) = math_inv3x3(crystallite_Fp(1:3,1:3,g,i,e))
            crystallite_Lp(1:3,1:3,g,i,e) = crystallite_subLp0(1:3,1:3,g,i,e)                           ! ...plastic velocity grad
            constitutive_state(g,i,e)%p = constitutive_subState0(g,i,e)%p                               ! ...microstructure
            crystallite_Tstar_v(1:6,g,i,e) = crystallite_subTstar0_v(1:6,g,i,e)                         ! ...2nd PK stress
                                                                                                        ! cant restore dotState here, since not yet calculated in first cutback after initialization
            !$OMP FLUSH(crystallite_invFp)
#ifndef _OPENMP
            if (debug_verbosity > 4 &
                .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
              write(6,'(a,f10.8)') '<< CRYST >> cutback step in crystallite_stressAndItsTangent with new crystallite_subStep: ',&
                                     crystallite_subStep(g,i,e)
              write(6,*)
            endif
#endif
          endif

          ! --- prepare for integration ---
          
          !$OMP FLUSH(crystallite_subStep)
          crystallite_todo(g,i,e) = crystallite_subStep(g,i,e) > subStepMinCryst                        ! still on track or already done (beyond repair)
          if (crystallite_todo(g,i,e)) then
            crystallite_subF(1:3,1:3,g,i,e) = crystallite_subF0(1:3,1:3,g,i,e) &
                                            + crystallite_subStep(g,i,e) &
                                              * (crystallite_partionedF(1:3,1:3,g,i,e) - crystallite_partionedF0(1:3,1:3,g,i,e))
            !$OMP FLUSH(crystallite_subF)
            crystallite_Fe(1:3,1:3,g,i,e) = math_mul33x33(crystallite_subF(1:3,1:3,g,i,e), crystallite_invFp(1:3,1:3,g,i,e))
            crystallite_subdt(g,i,e) = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
            crystallite_converged(g,i,e) = .false.                                                      ! start out non-converged
          endif
          
        enddo
      enddo
    enddo
  !$OMP END PARALLEL DO

  ! --- integrate ---
  
  if (any(crystallite_todo)) then
    select case(numerics_integrator(numerics_integrationMode))
      case(1)
        call crystallite_integrateStateFPI()
      case(2)
        call crystallite_integrateStateEuler()
      case(3)
        call crystallite_integrateStateAdaptiveEuler()
      case(4)
        call crystallite_integrateStateRK4()
      case(5)
        call crystallite_integrateStateRKCK45()
    endselect
  endif
  
  NiterationCrystallite = NiterationCrystallite + 1
      
enddo                                                                                                   ! cutback loop


! --+>> CHECK FOR NON-CONVERGED CRYSTALLITES <<+--

!$OMP PARALLEL DO PRIVATE(myNgrains,invFp,Fe_guess,Tstar)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
    myNgrains = homogenization_Ngrains(mesh_element(3,e))
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
      do g = 1,myNgrains
        if (.not. crystallite_converged(g,i,e)) then                                                    ! respond fully elastically (might be not required due to becoming terminally ill anyway)
          invFp = math_inv3x3(crystallite_partionedFp0(1:3,1:3,g,i,e))
          Fe_guess = math_mul33x33(crystallite_partionedF(1:3,1:3,g,i,e), invFp)
          Tstar = math_Mandel6to33( math_mul66x6( 0.5_pReal*constitutive_homogenizedC(g,i,e), &
                                                  math_Mandel33to6( math_mul33x33(transpose(Fe_guess),Fe_guess) - math_I3 ) ) )
          crystallite_P(1:3,1:3,g,i,e) = math_mul33x33(Fe_guess,math_mul33x33(Tstar,transpose(invFp)))
        endif
#ifndef _OPENMP
        if (debug_verbosity > 4 &
            .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
          write (6,'(a,i5,x,i2,x,i3)') '<< CRYST >> central solution of cryst_StressAndTangent at el ip g ',e,i,g
          write (6,*) 
          write (6,'(a,/,3(12(x),3(f12.4,x)/))') '<< CRYST >> P / MPa', math_transpose3x3(crystallite_P(1:3,1:3,g,i,e)) / 1e6
          write (6,'(a,/,3(12(x),3(f14.9,x)/))') '<< CRYST >> Fp', math_transpose3x3(crystallite_Fp(1:3,1:3,g,i,e))
          write (6,'(a,/,3(12(x),3(f14.9,x)/))') '<< CRYST >> Lp', math_transpose3x3(crystallite_Lp(1:3,1:3,g,i,e))
          write (6,*) 
        endif
#endif
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO


! --+>> STIFFNESS CALCULATION <<+--

if(updateJaco) then                                                                                     ! Jacobian required
  numerics_integrationMode = 2_pInt
  
  ! --- BACKUP ---
  
  !$OMP PARALLEL DO PRIVATE(myNgrains,mySizeState,mySizeDotState)
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                  ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          mySizeState = constitutive_sizeState(g,i,e)                                                   ! number of state variables for this grain
          mySizeDotState = constitutive_sizeDotState(g,i,e)                                             ! number of dotStates for this grain
          constitutive_state_backup(g,i,e)%p(1:mySizeState) = &
                 constitutive_state(g,i,e)%p(1:mySizeState)                                             ! remember unperturbed, converged state, ...
          constitutive_dotState_backup(g,i,e)%p(1:mySizeDotState) = &
                 constitutive_dotState(g,i,e)%p(1:mySizeDotState)                                       ! ... dotStates, ...
    enddo; enddo; enddo
  !$OMP END PARALLEL DO
  Temperature_backup = crystallite_Temperature                                                          ! ... Temperature, ...
  F_backup = crystallite_subF                                                                           ! ... and kinematics
  Fp_backup = crystallite_Fp
  InvFp_backup = crystallite_invFp
  Fe_backup = crystallite_Fe
  Lp_backup = crystallite_Lp
  Tstar_v_backup = crystallite_Tstar_v
  P_backup = crystallite_P
  convergenceFlag_backup = crystallite_converged

  
  ! --- CALCULATE STATE AND STRESS FOR PERTURBATION ---
  
  dPdF_perturbation1 = crystallite_dPdF0                                                                ! initialize stiffness with known good values from last increment
  dPdF_perturbation2 = crystallite_dPdF0                                                                ! initialize stiffness with known good values from last increment
  do perturbation = 1,2                                                                                 ! forward and backward perturbation
    if (iand(pert_method,perturbation) > 0) then                                                        ! mask for desired direction
      myPert = -pert_Fg * (-1.0_pReal)**perturbation                                                    ! set perturbation step
      do k = 1,3; do l = 1,3                                                                            ! ...alter individual components
        if (debug_verbosity> 5) then
          !$OMP CRITICAL (write2out)
            write(6,'(a,2(x,i1),x,a)') '<< CRYST >> [[[[[[ Stiffness perturbation',k,l,']]]]]]'
            write(6,*)
          !$OMP END CRITICAL (write2out)
        endif
        crystallite_subF(k,l,:,:,:) = crystallite_subF(k,l,:,:,:) + myPert                              ! perturb either forward or backward
        
        crystallite_todo = crystallite_requested .and. crystallite_converged
        where (crystallite_todo) crystallite_converged = .false.                                        ! start out non-converged

        select case(numerics_integrator(numerics_integrationMode))
          case(1)
            call crystallite_integrateStateFPI()
          case(2)
            call crystallite_integrateStateEuler()
          case(3)
            call crystallite_integrateStateAdaptiveEuler()
          case(4)
            call crystallite_integrateStateRK4()
          case(5)
            call crystallite_integrateStateRKCK45()
        end select
        
        !OMP PARALLEL DO PRIVATE(myNgrains)
          do e = FEsolving_execElem(1),FEsolving_execElem(2)
            myNgrains = homogenization_Ngrains(mesh_element(3,e))
            do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
              do g = 1,myNgrains
                if (crystallite_requested(g,i,e) .and. crystallite_converged(g,i,e)) then               ! converged state warrants stiffness update
                  select case(perturbation)
                    case(1)
                      dPdF_perturbation1(1:3,1:3,k,l,g,i,e) = (crystallite_P(1:3,1:3,g,i,e) - P_backup(1:3,1:3,g,i,e)) / myPert ! tangent dP_ij/dFg_kl
                    case(2)
                      dPdF_perturbation2(1:3,1:3,k,l,g,i,e) = (crystallite_P(1:3,1:3,g,i,e) - P_backup(1:3,1:3,g,i,e)) / myPert ! tangent dP_ij/dFg_kl
                  end select
                endif
          enddo; enddo; enddo
        !OMP END PARALLEL DO
        

        ! --- RESTORE ---
        
        !$OMP PARALLEL DO PRIVATE(myNgrains,mySizeState,mySizeDotState)
          do e = FEsolving_execElem(1),FEsolving_execElem(2)
            myNgrains = homogenization_Ngrains(mesh_element(3,e))
            do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
              do g = 1,myNgrains
                mySizeState = constitutive_sizeState(g,i,e)
                mySizeDotState = constitutive_sizeDotState(g,i,e)
                constitutive_state(g,i,e)%p(1:mySizeState) = constitutive_state_backup(g,i,e)%p(1:mySizeState)
                constitutive_dotState(g,i,e)%p(1:mySizeDotState) = constitutive_dotState_backup(g,i,e)%p(1:mySizeDotState)
          enddo; enddo; enddo
        !OMP END PARALLEL DO
        crystallite_Temperature = Temperature_backup
        crystallite_subF = F_backup
        crystallite_Fp = Fp_backup 
        crystallite_invFp = InvFp_backup
        crystallite_Fe = Fe_backup
        crystallite_Lp = Lp_backup
        crystallite_Tstar_v = Tstar_v_backup
        crystallite_P = P_backup
        crystallite_converged = convergenceFlag_backup

      enddo; enddo                                                                                      ! k,l loop
    endif
  enddo                                                                                                 ! perturbation direction


  ! --- STIFFNESS ACCORDING TO PERTURBATION METHOD AND CONVERGENCE ---

  !$OMP PARALLEL DO PRIVATE(myNgrains)
    do e = FEsolving_execElem(1),FEsolving_execElem(2)
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
        do g = 1,myNgrains
          if (crystallite_requested(g,i,e) .and. crystallite_converged(g,i,e)) then                     ! central solution converged
            select case(pert_method)
              case(1)
                crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = dPdF_perturbation1(1:3,1:3,1:3,1:3,g,i,e)
              case(2)
                crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = dPdF_perturbation2(1:3,1:3,1:3,1:3,g,i,e)
              case(3)
                crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = 0.5_pReal* (  dPdF_perturbation1(1:3,1:3,1:3,1:3,g,i,e) &
                                                                      + dPdF_perturbation2(1:3,1:3,1:3,1:3,g,i,e))
            end select
          elseif (crystallite_requested(g,i,e) .and. .not. crystallite_converged(g,i,e)) then           ! central solution did not converge
            crystallite_dPdF(1:3,1:3,1:3,1:3,g,i,e) = crystallite_fallbackdPdF(1:3,1:3,1:3,1:3,g,i,e)   ! use (elastic) fallback
          endif
    enddo; enddo; enddo
  !OMP END PARALLEL DO
  
endif                                                                                                   ! jacobian calculation
 
endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 4h order explicit Runge Kutta method 
!********************************************************************
subroutine crystallite_integrateStateRK4(gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use numerics, only:     numerics_integrationMode
use debug, only:        debug_verbosity, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_selectiveDebugger, &
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
                        constitutive_dotTemperature, &
                        constitutive_microstructure

implicit none

real(pReal), dimension(4), parameter ::       timeStepFraction = (/0.5_pReal, 0.5_pReal, 1.0_pReal, 1.0_pReal/) ! weight of slope used for Runge Kutta integration
real(pReal), dimension(4), parameter ::       weight = (/1.0_pReal, 2.0_pReal, 2.0_pReal, 1.0_pReal/)           ! factor giving the fraction of the original timestep used for Runge Kutta Integration

!*** input variables ***!
integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                              ii, &                     ! integration point index
                                              gg                        ! grain index

!*** output variables ***!

!*** local variables ***!
integer(pInt)                                 e, &                      ! element index in element loop
                                              i, &                      ! integration point index in ip loop
                                              g, &                      ! grain index in grain loop
                                              n, &
                                              mySizeDotState
integer(pInt), dimension(2) ::                eIter                     ! bounds for element iteration
integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                  ! bounds for ip iteration
                                              gIter                     ! bounds for grain iteration
real(pReal), dimension(homogenization_maxNgrains,mesh_maxNips,mesh_NcpElems) :: &
                                              RK4dotTemperature         ! evolution of Temperature of each grain for Runge Kutta integration
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
    gIter(1:2,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- RESET DOTSTATE ---

!$OMP PARALLEL PRIVATE(mySizeDotState)

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
 enddo; enddo; enddo
!$OMP ENDDO


! --- FIRST RUNGE KUTTA STEP ---

RK4dotTemperature = 0.0_pReal                                                                             ! initialize Runge-Kutta dotTemperature
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    constitutive_RK4dotState(g,i,e)%p = 0.0_pReal                                                         ! initialize Runge-Kutta dotState
    if (crystallite_todo(g,i,e)) then
      call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
    endif
  enddo; enddo; enddo
!$OMP ENDDO
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO


! --- SECOND TO FOURTH RUNGE KUTTA STEP PLUS FINAL INTEGRATION ---

do n = 1,4

  ! --- state update ---

  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        mySizeDotState = constitutive_sizeDotState(g,i,e)
        if (n < 4) then
          constitutive_RK4dotState(g,i,e)%p = constitutive_RK4dotState(g,i,e)%p + weight(n)*constitutive_dotState(g,i,e)%p
          RK4dotTemperature(g,i,e) = RK4dotTemperature(g,i,e) + weight(n)*crystallite_dotTemperature(g,i,e)
        elseif (n == 4) then
          constitutive_dotState(g,i,e)%p = (constitutive_RK4dotState(g,i,e)%p + weight(n)*constitutive_dotState(g,i,e)%p) /6.0_pReal  ! use weighted RKdotState for final integration
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
      endif
    enddo; enddo; enddo
  !$OMP ENDDO      

  
  ! --- update dependent states ---

  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)         ! update dependent state variables to be consistent with basic states
      endif
    enddo; enddo; enddo
  !$OMP ENDDO


  ! --- stress integration ---

  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        if (crystallite_integrateStress(g,i,e,timeStepFraction(n))) then                                  ! fraction of original times step
          if (n == 4) then                                                                                ! final integration step
#ifndef _OPENMP
            if (debug_verbosity > 5 &
                .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
              mySizeDotState = constitutive_sizeDotState(g,i,e)
              write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
              write(6,*)
              write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
              write(6,*)
              write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
              write(6,*)
            endif
#endif
            crystallite_converged(g,i,e) = .true.                                                         ! ... converged per definition
            crystallite_todo(g,i,e) = .false.                                                             ! ... integration done
            if (debug_verbosity > 0) then
              !$OMP CRITICAL (distributionState)
                debug_StateLoopDistribution(n,numerics_integrationMode) = &
                  debug_StateLoopDistribution(n,numerics_integrationMode) + 1
              !$OMP END CRITICAL (distributionState)
            endif
          endif
        else                                                                                              ! broken stress integration
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL (checkTodo)
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMP END CRITICAL (checkTodo)
          else                                                                                            ! if broken local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
          endif
        endif
      endif
      constitutive_dotState(g,i,e)%p = 0.0_pReal                                                          ! reset dotState to zero
    enddo; enddo; enddo
  !$OMP ENDDO      

  
  ! --- dot state and RK dot state---

  if (n < 4) then
    !$OMP DO
      do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
        if (crystallite_todo(g,i,e)) then
          call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                            crystallite_Temperature(g,i,e), timeStepFraction(n)*crystallite_subdt(g,i,e), & ! fraction of original timestep
                                            crystallite_orientation, g,i,e)
          crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                          crystallite_Temperature(g,i,e),g,i,e)
        endif
      enddo; enddo; enddo
    !$OMP ENDDO
    !$OMP DO
      do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                ! iterate over elements, ips and grains
        if (crystallite_todo(g,i,e)) then
          if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                    ! NaN occured in dotState
               .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then         ! NaN occured in dotTemperature
            if (.not. crystallite_localConstitution(g,i,e)) then                                          ! if broken non-local...
              !$OMP CRITICAL (checkTodo)
                crystallite_todo = crystallite_todo .and. crystallite_localConstitution                   ! ...all non-locals skipped
              !$OMP END CRITICAL (checkTodo)
            else                                                                                          ! if broken local...
              crystallite_todo(g,i,e) = .false.                                                           ! ... skip this one next time
            endif
          endif
        endif
      enddo; enddo; enddo
    !$OMP ENDDO
  endif
  
enddo

!$OMP END PARALLEL


! --- CHECK CONVERGENCE ---

crystallite_todo = .false.                                                                                ! done with integration
if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
  if (any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                    ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 5th order Runge-Kutta Cash-Karp method with adaptive step size
! (use 5th order solution to advance = "local extrapolation")
!********************************************************************
subroutine crystallite_integrateStateRKCK45(gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debug_verbosity, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_selectiveDebugger, &
                        debug_StateLoopDistribution
use numerics, only:     rTol_crystalliteState, &
                        rTol_crystalliteTemperature, &
                        subStepSizeCryst, &
                        stepIncreaseCryst, &
                        numerics_integrationMode
use FEsolving, only:    FEsolving_execElem, & 
                        FEsolving_execIP, &
                        theInc
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
                        constitutive_dotTemperature, &
                        constitutive_microstructure

implicit none


!*** input variables ***!
integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                              ii, &                     ! integration point index
                                              gg                        ! grain index

!*** output variables ***!

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
    gIter(1:2,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- RESET DOTSTATE ---

!$OMP PARALLEL PRIVATE(mySizeDotState)

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
  enddo; enddo; enddo
!$OMP ENDDO


! --- FIRST RUNGE KUTTA STEP ---
#ifndef _OPENMP
if (debug_verbosity > 5) then
  write(6,'(a,x,i1)') '<< CRYST >> RUNGE KUTTA STEP',1
endif
#endif
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
    endif
  enddo; enddo; enddo
!$OMP ENDDO
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      endif
    endif
 enddo; enddo; enddo
!$OMP ENDDO


! --- SECOND TO SIXTH RUNGE KUTTA STEP ---

do n = 1,5

  ! --- state update ---
  
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

  
  ! --- update dependent states ---

  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)         ! update dependent state variables to be consistent with basic states
      endif
      constitutive_dotState(g,i,e)%p = 0.0_pReal                                                          ! reset dotState to zero
    enddo; enddo; enddo
  !$OMP ENDDO


  ! --- stress integration ---
  
  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        if (.not. crystallite_integrateStress(g,i,e,c(n))) then                                           ! fraction of original time step
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL (checkTodo)
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMP END CRITICAL (checkTodo)
          else                                                                                            ! if broken local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
          endif
        endif
      endif
    enddo; enddo; enddo
  !$OMP ENDDO      


  ! --- dot state and RK dot state---
#ifndef _OPENMP
  if (debug_verbosity > 5) then
    write(6,'(a,x,i1)') '<< CRYST >> RUNGE KUTTA STEP',n+1
  endif
#endif
  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                          crystallite_Temperature(g,i,e), c(n)*crystallite_subdt(g,i,e), & ! fraction of original timestep
                                          crystallite_orientation, g,i,e)
        crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                        crystallite_Temperature(g,i,e),g,i,e)
      endif
    enddo; enddo; enddo
  !$OMP ENDDO
  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
             .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL (checkTodo)
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMP END CRITICAL (checkTodo)
          else                                                                                            ! if broken local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
          endif
        endif
      endif
    enddo; enddo; enddo
  !$OMP ENDDO

enddo  


! --- STATE UPDATE WITH ERROR ESTIMATE FOR STATE AND TEMPERATURE ---

relStateResiduum = 0.0_pReal
relTemperatureResiduum = 0.0_pReal
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
      
      stateResiduum(1:mySizeDotState,g,i,e) = (   db(1) * constitutive_RKCK45dotState(1,g,i,e)%p(1:mySizeDotState) &
                                                + db(2) * constitutive_RKCK45dotState(2,g,i,e)%p(1:mySizeDotState) &
                                                + db(3) * constitutive_RKCK45dotState(3,g,i,e)%p(1:mySizeDotState) &
                                                + db(4) * constitutive_RKCK45dotState(4,g,i,e)%p(1:mySizeDotState) &
                                                + db(5) * constitutive_RKCK45dotState(5,g,i,e)%p(1:mySizeDotState) &
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
      forall (s = 1:mySizeDotState, abs(constitutive_state(g,i,e)%p(s)) > 0.0_pReal) &
        relStateResiduum(s,g,i,e) = stateResiduum(s,g,i,e) / constitutive_state(g,i,e)%p(s)
      if (crystallite_Temperature(g,i,e) > 0) &
        relTemperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) / crystallite_Temperature(g,i,e)
      
      !$OMP FLUSH(relStateResiduum,relTemperatureResiduum)
      
      crystallite_todo(g,i,e) = &
          ( all(      abs(relStateResiduum(:,g,i,e)) < rTol_crystalliteState &
                 .or. abs(stateResiduum(1:mySizeDotState,g,i,e)) < constitutive_aTolState(g,i,e)%p(1:mySizeDotState) ) &
           .and. abs(relTemperatureResiduum(g,i,e)) < rTol_crystalliteTemperature )
           
#ifndef _OPENMP
      if (debug_verbosity > 5 &
          .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
        write(6,'(a,i5,x,i3,x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
        write(6,*)
        write(6,'(a,/,(12(x),12(f12.1,x)))') '<< CRYST >> absolute residuum tolerance', &
                                        stateResiduum(1:mySizeDotState,g,i,e) / constitutive_aTolState(g,i,e)%p(1:mySizeDotState)
        write(6,*)
        write(6,'(a,/,(12(x),12(f12.1,x)))') '<< CRYST >> relative residuum tolerance', &
                                              relStateResiduum(1:mySizeDotState,g,i,e) / rTol_crystalliteState
        write(6,*)
        write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
        write(6,*)
        write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
        write(6,*)
      endif
#endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO
      

! --- UPDATE DEPENDENT STATES IF RESIDUUM BELOW TOLERANCE ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMP ENDDO


! --- FINAL STRESS INTEGRATION STEP IF RESIDUUM BELOW TOLERANCE ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if (crystallite_integrateStress(g,i,e)) then
        crystallite_converged(g,i,e) = .true.                                                           ! ... converged per definitionem
        crystallite_todo(g,i,e) = .false.                                                               ! ... integration done
        if (debug_verbosity > 0) then
          !$OMP CRITICAL (distributionState)
            debug_StateLoopDistribution(6,numerics_integrationMode) = debug_StateLoopDistribution(6,numerics_integrationMode) + 1
          !$OMP END CRITICAL (distributionState)
        endif
      else
        if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        endif        
      endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO

!$OMP END PARALLEL

! --- nonlocal convergence check ---

#ifndef _OPENMP  
  if (debug_verbosity > 5) then
    write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), ' grains converged'
    write(6,*)
  endif
#endif
if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
  if ( any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                   ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
  
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 1nd order Euler method with adaptive step size
!********************************************************************
subroutine crystallite_integrateStateAdaptiveEuler(gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debug_verbosity, &
                        debug_selectiveDebugger, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_StateLoopDistribution
use numerics, only:     rTol_crystalliteState, &
                        rTol_crystalliteTemperature, &
                        subStepSizeCryst, &
                        stepIncreaseCryst, &
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

!*** output variables ***!

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
    gIter(1:2,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- RESET DOTSTATE ---

!$OMP PARALLEL PRIVATE(mySizeDotState)

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
 enddo; enddo; enddo
!$OMP ENDDO


! --- DOT STATE AND TEMPERATURE (EULER INTEGRATION) ---

stateResiduum = 0.0_pReal
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then  
      call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
    endif
 enddo; enddo; enddo
!$OMP ENDDO
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then  
      if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
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
      constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                                    + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
      crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                     + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
    endif
  enddo; enddo; enddo
!$OMP ENDDO      


! --- UPDATE DEPENDENT STATES (EULER INTEGRATION) ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)           ! update dependent state variables to be consistent with basic states
    endif
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
 enddo; enddo; enddo
!$OMP ENDDO


! --- STRESS INTEGRATION (EULER INTEGRATION) ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if (.not. crystallite_integrateStress(g,i,e)) then
        if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        else                                                                                            ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO      


! --- DOT STATE AND TEMPERATURE (HEUN METHOD) ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
    endif
  enddo; enddo; enddo
!$OMP ENDDO
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if ( any(constitutive_dotState(g,i,e)%p /= constitutive_dotState(g,i,e)%p) &                      ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e) /= crystallite_dotTemperature(g,i,e) ) then           ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        else                                                                                            ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
        endif
      endif      
    endif
  enddo; enddo; enddo
!$OMP ENDDO


! --- ERROR ESTIMATE FOR STATE AND TEMPERATURE (HEUN METHOD) ---

relStateResiduum = 0.0_pReal
relTemperatureResiduum = 0.0_pReal
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

      forall (s = 1:mySizeDotState, abs(constitutive_state(g,i,e)%p(s)) > 0.0_pReal) &
        relStateResiduum(s,g,i,e) = stateResiduum(s,g,i,e) / constitutive_state(g,i,e)%p(s)
      if (crystallite_Temperature(g,i,e) > 0) &
        relTemperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) / crystallite_Temperature(g,i,e)
      !$OMP FLUSH(relStateResiduum,relTemperatureResiduum)

#ifndef _OPENMP        
      if (debug_verbosity > 5 &
          .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
        write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
        write(6,*)
        write(6,'(a,/,(12(x),12(f12.1,x)))') '<< CRYST >> absolute residuum tolerance', &
                                        stateResiduum(1:mySizeDotState,g,i,e) / constitutive_aTolState(g,i,e)%p(1:mySizeDotState)
        write(6,*)
        write(6,'(a,/,(12(x),12(f12.1,x)))') '<< CRYST >> relative residuum tolerance', &
                                              relStateResiduum(1:mySizeDotState,g,i,e) / rTol_crystalliteState
        write(6,*)
        write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState) &
                                                    - 2.0_pReal * stateResiduum(1:mySizeDotState,g,i,e) / crystallite_subdt(g,i,e)  ! calculate former dotstate from higher order solution and state residuum
        write(6,*)
        write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
        write(6,*)
      endif
#endif            
      
      ! --- converged ? ---

      if ( all(     abs(relStateResiduum(:,g,i,e)) < rTol_crystalliteState &
               .or. abs(stateResiduum(1:mySizeDotState,g,i,e)) < constitutive_aTolState(g,i,e)%p(1:mySizeDotState)) &
          .and. abs(relTemperatureResiduum(g,i,e)) < rTol_crystalliteTemperature ) then        
        crystallite_converged(g,i,e) = .true.                                                             ! ... converged per definitionem
        crystallite_todo(g,i,e) = .false.                                                                 ! ... integration done
        if (debug_verbosity > 0) then
          !$OMP CRITICAL (distributionState)
            debug_StateLoopDistribution(2,numerics_integrationMode) = debug_StateLoopDistribution(2,numerics_integrationMode) + 1
          !$OMP END CRITICAL (distributionState)
        endif
      endif

    endif
  enddo; enddo; enddo
!$OMP ENDDO

!$OMP END PARALLEL

! --- NONLOCAL CONVERGENCE CHECK ---

#ifndef _OPENMP  
  if (debug_verbosity > 5) then
    write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), ' grains converged'
    write(6,*)
  endif
#endif
if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
  if ( any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                   ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 1st order explicit Euler method 
!********************************************************************
subroutine crystallite_integrateStateEuler(gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use numerics, only:     numerics_integrationMode
use debug, only:        debug_verbosity, &
                        debug_selectiveDebugger, &
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

!*** output variables ***!

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
    gIter(1:2,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- RESET DOTSTATE ---

!$OMP PARALLEL PRIVATE(mySizeDotState)

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
 enddo; enddo; enddo
!$OMP ENDDO


! --- DOT STATE AND TEMPERATURE ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
    endif
  enddo; enddo; enddo
!$OMP ENDDO
!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO


! --- UPDATE STATE AND TEMPERATURE ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      mySizeDotState = constitutive_sizeDotState(g,i,e)
      constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                                    + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
      crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                        + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
#ifndef _OPENMP
      if (debug_verbosity > 5 &
          .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
        write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateState at el ip g ',e,i,g
        write(6,*)
        write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
        write(6,*)
        write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
        write(6,*)
      endif
#endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO

  
! --- UPDATE DEPENDENT STATES ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMP ENDDO


! --- STRESS INTEGRATION ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      if (crystallite_integrateStress(g,i,e)) then
        crystallite_converged(g,i,e) = .true.
        if (debug_verbosity > 0) then
          !$OMP CRITICAL (distributionState)
            debug_StateLoopDistribution(1,numerics_integrationMode) = debug_StateLoopDistribution(1,numerics_integrationMode) + 1
          !$OMP END CRITICAL (distributionState)
        endif
      else                                                                                                ! broken stress integration
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL (checkTodo)
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO      

!$OMP END PARALLEL


! --- CHECK NON-LOCAL CONVERGENCE ---

crystallite_todo = .false.                                                                                ! done with integration
if (.not. singleRun) then                                                                                 ! if not requesting Integration of just a single IP   
  if (any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                    ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! adaptive 1st order explicit Euler method 
! using Fixed Point Iteration to adapt the stepsize  
!********************************************************************
subroutine crystallite_integrateStateFPI(gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debug_verbosity, &
                        debug_selectiveDebugger, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_StateLoopDistribution
use numerics, only:     nState, &
                        numerics_integrationMode
use FEsolving, only:    FEsolving_execElem, & 
                        FEsolving_execIP
use mesh, only:         mesh_element, &
                        mesh_NcpElems
use material, only:     homogenization_Ngrains
use constitutive, only: constitutive_sizeDotState, &
                        constitutive_state, &
                        constitutive_dotState, &
                        constitutive_collectDotState, &
                        constitutive_dotTemperature, &
                        constitutive_microstructure, &
                        constitutive_previousDotState, &
                        constitutive_previousDotState2

implicit none

!*** input variables ***!
integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                              ii, &                     ! integration point index
                                              gg                        ! grain index

!*** output variables ***!

!*** local variables ***!
integer(pInt)                                 NiterationState, &        ! number of iterations in state loop
                                              e, &                      ! element index in element loop
                                              i, &                      ! integration point index in ip loop
                                              g                         ! grain index in grain loop
integer(pInt), dimension(2) ::                eIter                     ! bounds for element iteration
integer(pInt), dimension(2,mesh_NcpElems) ::  iIter, &                  ! bounds for ip iteration
                                              gIter                     ! bounds for grain iteration
real(pReal)                                   dot_prod12, &
                                              dot_prod22
logical                                       singleRun, &              ! flag indicating computation for single (g,i,e) triple
                                              stateConverged, &         ! flag indicating convergence of state integration
                                              temperatureConverged, &   ! flag indicating convergence of temperature integration
                                              stateUpdateDone, &        ! flag indicating successfull state update
                                              temperatureUpdateDone     ! flag indicating successfull temperature update


if (present(ee) .and. present(ii) .and. present(gg)) then
  eIter = ee
  iIter(1:2,ee) = ii
  gIter(1:2,ee) = gg
  singleRun = .true.
else
  eIter = FEsolving_execElem(1:2)
  do e = eIter(1),eIter(2)
    iIter(1:2,e) = FEsolving_execIP(1:2,e)
    gIter(1:2,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --+>> PREGUESS FOR STATE <<+--

! --- RESET DOTSTATE ---

!$OMP PARALLEL

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
    constitutive_previousDotState(g,i,e)%p = 0.0_pReal
    constitutive_previousDotState2(g,i,e)%p = 0.0_pReal
 enddo; enddo; enddo
!$OMP ENDDO


! --- DOT STATES ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
    endif
  enddo; enddo; enddo
!$OMP ENDDO


! --- STATE & TEMPERATURE UPDATE ---

!$OMP SINGLE
  crystallite_statedamper = 1.0_pReal
!$OMP END SINGLE
!$OMP DO PRIVATE(stateUpdateDone,temperatureUpdateDone,stateConverged,temperatureConverged)
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call crystallite_updateState(stateUpdateDone, stateConverged, g,i,e)                                ! update state
      call crystallite_updateTemperature(temperatureUpdateDone, temperatureConverged, g,i,e)              ! update temperature
      crystallite_todo(g,i,e) = stateUpdateDone .and. temperatureUpdateDone
      if ( (.not. stateUpdateDone .or. .not. temperatureUpdateDone) &
          .and. .not. crystallite_localConstitution(g,i,e) ) then                                         ! if updateState or updateTemperature signals broken non-local... 
        !$OMP CRITICAL (checkTodo) 
          crystallite_todo = crystallite_todo .and. crystallite_localConstitution                         ! ...all non-locals skipped
        !$OMP END CRITICAL (checkTodo)
      endif
    endif
  enddo; enddo; enddo
!$OMP ENDDO


! --- UPDATE DEPENDENT STATES ---

!$OMP DO
  do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                    ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)           ! update dependent state variables to be consistent with basic states
    endif
    constitutive_previousDotState2(g,i,e)%p = constitutive_previousDotState(g,i,e)%p                      ! age previous dotState
    constitutive_previousDotState(g,i,e)%p = constitutive_dotState(g,i,e)%p                               ! age previous dotState
    constitutive_dotState(g,i,e)%p = 0.0_pReal                                                            ! reset dotState to zero
  enddo; enddo; enddo
!$OMP ENDDO
!$OMP END PARALLEL


! --+>> STATE LOOP <<+--

NiterationState = 0_pInt
do while (any(crystallite_todo) .and. NiterationState < nState )                                          ! convergence loop for crystallite
  NiterationState = NiterationState + 1_pInt
  
  !$OMP PARALLEL
  

  ! --- STRESS INTEGRATION ---
  
  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        if (.not. crystallite_integrateStress(g,i,e)) then                                                ! if broken ...
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! ... and non-local... 
            !$OMP CRITICAL (checkTodo) 
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ... then all non-locals skipped
            !$OMP END CRITICAL (checkTodo)
          else                                                                                            ! ... and local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... then skip only me
          endif
        endif
      endif
    enddo; enddo; enddo
  !$OMP ENDDO

#ifndef _OPENMP
  if (debug_verbosity > 5) then
    write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo after stress integration'
  endif
#endif


  ! --- DOT STATES ---
  
  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_collectDotState(crystallite_Tstar_v(1:6,g,i,e), crystallite_Fe, crystallite_Fp, &
                                        crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), crystallite_orientation, g,i,e)
      endif
  enddo; enddo; enddo
  !$OMP ENDDO


  ! --- STATE & TEMPERATURE UPDATE ---

  !$OMP SINGLE
    crystallite_statedamper = 1.0_pReal
  !$OMP END SINGLE
  !$OMP DO PRIVATE(dot_prod12,dot_prod22,stateUpdateDone,temperatureUpdateDone,stateConverged,temperatureConverged)
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then

        ! --- state damper ---
        
        dot_prod12 = dot_product( constitutive_dotState(g,i,e)%p - constitutive_previousDotState(g,i,e)%p, &
                              constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p )
        dot_prod22 = dot_product( constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p, &
                                  constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p )
        if (      dot_prod22 > 0.0_pReal &
            .and. (     dot_prod12 < 0.0_pReal &
                   .or. dot_product(constitutive_dotState(g,i,e)%p, constitutive_previousDotState(g,i,e)%p) < 0.0_pReal) ) then
          crystallite_statedamper(g,i,e) = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)
          !$OMP FLUSH(crystallite_statedamper)
        endif
        
        ! --- updates ---
        
        call crystallite_updateState(stateUpdateDone, stateConverged, g,i,e)                              ! update state
        call crystallite_updateTemperature(temperatureUpdateDone, temperatureConverged, g,i,e)            ! update temperature
        crystallite_todo(g,i,e) = stateUpdateDone .and. temperatureUpdateDone
        crystallite_converged(g,i,e) = stateConverged .and. temperatureConverged
        if ( (.not. stateUpdateDone .or. .not. temperatureUpdateDone) &
            .and. .not. crystallite_localConstitution(g,i,e) ) then                                       ! if updateState or updateTemperature signals broken non-local... 
          !$OMP CRITICAL (checkTodo) 
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMP END CRITICAL (checkTodo)
        elseif (stateConverged .and. temperatureConverged) then                                           ! check (private) logicals "stateConverged" and "temperatureConverged" instead of (shared) "crystallite_converged", so no need to flush the "crystallite_converged" array
          if (debug_verbosity > 0) then
            !$OMP CRITICAL (distributionState)
              debug_StateLoopDistribution(NiterationState,numerics_integrationMode) = &
                debug_StateLoopDistribution(NiterationState,numerics_integrationMode) + 1
            !$OMP END CRITICAL (distributionState)
          endif
        endif
      endif
    enddo; enddo; enddo
  !$OMP ENDDO


  ! --- UPDATE DEPENDENT STATES ---

  !$OMP DO
    do e = eIter(1),eIter(2); do i = iIter(1,e),iIter(2,e); do g = gIter(1,e),gIter(2,e)                  ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Fe, g, i, e)         ! update dependent state variables to be consistent with basic states
      endif
      constitutive_previousDotState2(g,i,e)%p = constitutive_previousDotState(g,i,e)%p                    ! age previous dotState
      constitutive_previousDotState(g,i,e)%p = constitutive_dotState(g,i,e)%p                             ! age previous dotState
      constitutive_dotState(g,i,e)%p = 0.0_pReal                                                          ! reset dotState to zero
   enddo; enddo; enddo
  !$OMP ENDDO
  
  !$OMP END PARALLEL

#ifndef _OPENMP  
  if (debug_verbosity > 5) then
    write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_converged(:,:,:)), &
                              ' grains converged after state integration no. ', NiterationState
    write(6,*)
  endif
#endif

  
  ! --- CONVERGENCE CHECK ---

  if (.not. singleRun) then                                                                               ! if not requesting Integration of just a single IP   
    if (any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                  ! any non-local not yet converged (or broken)...
      crystallite_converged = crystallite_converged .and. crystallite_localConstitution                   ! ...restart all non-local as not converged
    endif
  endif
  crystallite_todo = crystallite_todo .and. .not. crystallite_converged                                   ! skip all converged
  
#ifndef _OPENMP
  if (debug_verbosity > 5) then
    write(6,'(a,i8,a)') '<< CRYST >> ', count(crystallite_converged(:,:,:)),' grains converged after non-local check'
    write(6,'(a,i8,a,i2)') '<< CRYST >> ', count(crystallite_todo(:,:,:)),' grains todo after state integration no. ',&
                                NiterationState
    write(6,*)
  endif
#endif
  
enddo                                                                                           ! crystallite convergence loop  

endsubroutine



!********************************************************************
! update the internal state of the constitutive law
! and tell whether state has converged
!********************************************************************
subroutine crystallite_updateState(done, converged, g, i, e)

!*** variables and functions from other modules ***!
use prec, only:           pReal, &
                          pInt, &
                          pLongInt
use numerics, only:       rTol_crystalliteState, &
                          numerics_integrationMode
use constitutive, only:   constitutive_dotState, &
                          constitutive_previousDotState, &
                          constitutive_sizeDotState, &
                          constitutive_subState0, &
                          constitutive_state, &
                          constitutive_aTolState, &
                          constitutive_microstructure
use debug, only:          debug_verbosity, &
                          debug_g, &
                          debug_i, &
                          debug_e, &
                          debug_selectiveDebugger

!*** input variables ***!
integer(pInt), intent(in):: e, &                      ! element index
                            i, &                      ! integration point index
                            g                         ! grain index

!*** output variables ***!
logical, intent(out) ::     converged, &              ! flag indicating if state converged
                            done                      ! flag indicating if state was updated

!*** local variables ***!
real(pReal), dimension(constitutive_sizeDotState(g,i,e)) :: &
                            residuum, &              ! residuum from evolution of microstructure
                            dotState, &              ! local copy of dotState
                            state                    ! local copy of relevant part of state
integer(pInt)               mySize


!* start out as not done and not converged 
!* and make local copies of state and dotState (needed for parallelization, since no chance of flushing a pointer)

done = .false.
converged = .false.
mySize = constitutive_sizeDotState(g,i,e)
dotState(1:mySize) = constitutive_dotState(g,i,e)%p(1:mySize)
state(1:mySize) = constitutive_state(g,i,e)%p(1:mySize)


!* correct dotState, calculate residuum and check if it is valid

dotState(1:mySize) = constitutive_dotState(g,i,e)%p(1:mySize) * crystallite_statedamper(g,i,e) &
                   + constitutive_previousDotState(g,i,e)%p(1:mySize) * (1.0_pReal - crystallite_statedamper(g,i,e))
residuum = constitutive_state(g,i,e)%p(1:mySize) - constitutive_subState0(g,i,e)%p(1:mySize) &
                                                 - dotState(1:mySize) * crystallite_subdt(g,i,e)
if (any(residuum /= residuum)) then                                   ! if NaN occured then return without changing the state
#ifndef _OPENMP
  if (debug_verbosity > 4) then
    write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateState encountered NaN at el ip g ',e,i,g
  endif
#endif
  return
endif


!* update state and set convergence flag to true if residuum is below relative/absolute tolerance, otherwise set it to false

state(1:mySize) = state(1:mySize) - residuum
done = .true.
converged = all(     abs(residuum) < constitutive_aTolState(g,i,e)%p(1:mySize) &
                .or. abs(residuum) < rTol_crystalliteState * abs(state(1:mySize)) )

#ifndef _OPENMP
if (debug_verbosity > 5 .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
  if (converged) then
    write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateState converged at el ip g ',e,i,g
  else
    write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateState did not converge at el ip g ',e,i,g
  endif
  write(6,*)
  write(6,'(a,f6.1)') '<< CRYST >> crystallite_statedamper ',crystallite_statedamper(g,i,e)
  write(6,*)
  write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> dotState',dotState(1:mySize)
  write(6,*)
  write(6,'(a,/,(12(x),12(e12.5,x)))') '<< CRYST >> new state',state(1:mySize)
  write(6,*)
endif
#endif


!* sync global state and dotState with local copies

constitutive_dotState(g,i,e)%p(1:mySize) = dotState(1:mySize)
constitutive_state(g,i,e)%p(1:mySize) = state(1:mySize)

endsubroutine



!********************************************************************
! update the temperature of the grain
! and tell whether it has converged
!********************************************************************
 subroutine crystallite_updateTemperature(done, converged, g, i, e)
 
!*** variables and functions from other modules ***!
use prec, only:         pReal, &
                        pInt, &
                        pLongInt
use numerics, only:     rTol_crystalliteTemperature
use constitutive, only: constitutive_dotTemperature
use debug, only:        debug_verbosity
 
!*** input variables ***!
integer(pInt), intent(in):: e, &                      ! element index
                            i, &                      ! integration point index
                            g                         ! grain index
 
!*** output variables ***!
logical, intent(out) ::     converged, &              ! flag indicating if temperature converged
                            done                      ! flag indicating if temperature was updated

!*** local variables ***!
real(pReal)                 residuum                  ! residuum from evolution of temperature
 

!* start out as not done and not converged 

done = .false.
converged = .false.


!* correct dotState, calculate residuum and check if it is valid
!* (if NaN occured then return without changing the temperature)

residuum = crystallite_Temperature(g,i,e) - crystallite_subTemperature0(g,i,e) &
         - constitutive_dotTemperature(crystallite_Tstar_v(1:6,g,i,e),crystallite_Temperature(g,i,e),g,i,e) &
         * crystallite_subdt(g,i,e)
if (residuum /= residuum) then
#ifndef _OPENMP
  if (debug_verbosity > 4) then
    write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> updateTemperature encountered NaN at el ip g ',e,i,g
  endif
#endif
  return
endif
 

!* update temperature and set convergence flag to true if residuum is below relative tolerance (or zero Kelvin), otherwise set it to false

crystallite_Temperature(g,i,e) = crystallite_Temperature(g,i,e) - residuum
done = .true.
!$OMP FLUSH(crystallite_Temperature)
converged = (      crystallite_Temperature(g,i,e) == 0.0_pReal &
              .or. abs(residuum) < rTol_crystalliteTemperature * crystallite_Temperature(g,i,e))
 
endsubroutine



!***********************************************************************
!***     calculation of stress (P) with time integration             ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
function crystallite_integrateStress(&
     g,&          ! grain number
     i,&          ! integration point number
     e,&          ! element number
     fraction &
     )
     

!*** variables and functions from other modules ***!
use prec, only:         pReal, &
                        pInt, &
                        pLongInt
use numerics, only:     nStress, &
                        aTol_crystalliteStress, &
                        rTol_crystalliteStress, &
                        iJacoLpresiduum, &
                        relevantStrain, &
                        numerics_integrationMode
use debug, only:        debug_verbosity, &
                        debug_g, &
                        debug_i, &
                        debug_e, &
                        debug_selectiveDebugger, &
                        debug_cumLpCalls, &
                        debug_cumLpTicks, &
                        debug_StressLoopDistribution, &
                        debug_LeapfrogBreakDistribution
use constitutive, only: constitutive_homogenizedC, &
                        constitutive_LpAndItsTangent
use math, only:         math_mul33x33, &
                        math_mul66x6, &
                        math_mul99x99, &
                        math_transpose3x3, &
                        math_inv3x3, &
                        math_invert3x3, &
                        math_invert, &
                        math_det3x3, &
                        math_I3, &
                        math_identity2nd, &
                        math_Mandel66to3333, &
                        math_Mandel6to33, &
                        math_mandel33to6

implicit none

!*** input variables ***!
integer(pInt), intent(in)::         e, &                          ! element index
                                    i, &                          ! integration point index
                                    g                             ! grain index
real(pReal), optional, intent(in) :: fraction                      ! fraction of timestep

!*** output variables ***!
logical                             crystallite_integrateStress   ! flag indicating if integration suceeded

!*** local variables ***!
real(pReal), dimension(3,3)::       Fg_new, &                     ! deformation gradient at end of timestep
                                    Fp_current, &                 ! plastic deformation gradient at start of timestep
                                    Fp_new, &                     ! plastic deformation gradient at end of timestep
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
real(pReal), dimension(6)::         Tstar_v                       ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
real(pReal), dimension(9,9)::       dLpdT_constitutive, &         ! partial derivative of plastic velocity gradient calculated by constitutive law
                                    dTdLp, &                      ! partial derivative of 2nd Piola-Kirchhoff stress
                                    dRdLp, &                      ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                    invdRdLp                      ! inverse of dRdLp
real(pReal), dimension(3,3,3,3)::   C                             ! 4th rank elasticity tensor
real(pReal), dimension(6,6)::       C_66                          ! simplified 2nd rank elasticity tensor 
real(pReal)                         p_hydro, &                    ! volumetric part of 2nd Piola-Kirchhoff Stress
                                    det, &                        ! determinant
                                    leapfrog, &                   ! acceleration factor for Newton-Raphson scheme
                                    maxleap, &                    ! maximum acceleration factor
                                    dt                            ! time increment
logical                             error                         ! flag indicating an error
integer(pInt)                       NiterationStress, &           ! number of stress integrations
                                    dummy, &
                                    h, &
                                    j, &
                                    k, &
                                    l, &
                                    m, &
                                    n, &
                                    jacoCounter                   ! counter to check for Jacobian update
integer(pLongInt)                   tick, &
                                    tock, &
                                    tickrate, &
                                    maxticks
 
!* be pessimistic

crystallite_integrateStress = .false.
#ifndef _OPENMP
if (debug_verbosity > 5 .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
  write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> integrateStress at el ip g ',e,i,g
endif
#endif


!* only integrate over fraction of timestep?

if (present(fraction)) then
  dt = crystallite_subdt(g,i,e) * fraction
  Fg_new = crystallite_subF0(1:3,1:3,g,i,e) + (crystallite_subF(1:3,1:3,g,i,e) - crystallite_subF0(1:3,1:3,g,i,e)) * fraction
else     
  dt = crystallite_subdt(g,i,e)
  Fg_new = crystallite_subF(1:3,1:3,g,i,e)
endif

 
!* feed local variables

Fp_current =   crystallite_subFp0(1:3,1:3,g,i,e)
Tstar_v =      crystallite_Tstar_v(1:6,g,i,e)
Lpguess_old =  crystallite_Lp(1:3,1:3,g,i,e)                       ! consider present Lp good (i.e. worth remembering) ...
Lpguess =      crystallite_Lp(1:3,1:3,g,i,e)                       ! ... and take it as first guess


!* inversion of Fp_current...

invFp_current = math_inv3x3(Fp_current)                            
if (all(invFp_current == 0.0_pReal)) then                          ! ... failed?
#ifndef _OPENMP
  if (debug_verbosity > 4) then
    write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> integrateStress failed on invFp_current inversion at el ip g ',e,i,g
    if (debug_verbosity > 5 .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
      write(6,*)
      write(6,'(a,/,3(12(x),3(f12.7,x)/))') '<< CRYST >> invFp_new',math_transpose3x3(invFp_new(1:3,1:3))
    endif
  endif
#endif
  return
endif
A = math_mul33x33(transpose(invFp_current), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_current)))
 
 
!* get elasticity tensor

C_66 = constitutive_homogenizedC(g,i,e)
C = math_Mandel66to3333(C_66)


!* start LpLoop with no acceleration

NiterationStress = 0_pInt
leapfrog = 1.0_pReal
maxleap = 16.0_pReal
jacoCounter = 0_pInt

LpLoop: do
  NiterationStress = NiterationStress + 1
   
  
  !* too many loops required ?
  
  if (NiterationStress > nStress) then
#ifndef _OPENMP
    if (debug_verbosity > 4) then 
      write(6,'(a,i5,x,i2,x,i3)') '<< CRYST >> integrateStress reached loop limit at el ip g ',e,i,g
      write(6,*)
    endif
#endif
    return
  endif
   
  B = math_I3 - dt*Lpguess
  BT = math_transpose3x3(B)
  AB = math_mul33x33(A,B)
  BTA = math_mul33x33(BT,A)
   
  
  !* calculate 2nd Piola-Kirchhoff stress tensor
  
  Tstar_v = 0.5_pReal * math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB) - math_I3))
  p_hydro = sum(Tstar_v(1:3)) / 3.0_pReal
  forall(n=1:3) Tstar_v(n) = Tstar_v(n) - p_hydro                  ! get deviatoric stress tensor
   
  
  !* calculate plastic velocity gradient and its tangent according to constitutive law
  
  if (debug_verbosity > 0) then
    call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
  endif
  call constitutive_LpAndItsTangent(Lp_constitutive, dLpdT_constitutive, Tstar_v, crystallite_Temperature(g,i,e), g, i, e)
  if (debug_verbosity > 0) then
    call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
    !$OMP CRITICAL (debugTimingLpTangent)
      debug_cumLpCalls = debug_cumLpCalls + 1_pInt
      debug_cumLpTicks = debug_cumLpTicks + tock-tick
      !$OMP FLUSH (debug_cumLpTicks)
      if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks
    !$OMP END CRITICAL (debugTimingLpTangent)
  endif
   
#ifndef _OPENMP
  if (debug_verbosity > 5 .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger) &
      .and. numerics_integrationMode == 1_pInt) then
    write(6,'(a,i3)') '<< CRYST >> iteration ', NiterationStress
    write(6,*)
    write(6,'(a,/,3(12(x),3(e20.7,x)/))') '<< CRYST >> Lp_constitutive', math_transpose3x3(Lp_constitutive)
    write(6,'(a,/,3(12(x),3(e20.7,x)/))') '<< CRYST >> Lpguess', math_transpose3x3(Lpguess)
  endif
#endif


  !* update current residuum and check for convergence of loop
  
  residuum = Lpguess - Lp_constitutive
  if (.not.(any(residuum/=residuum)) .and. &                       ! exclude any NaN in residuum
       ( maxval(abs(residuum)) < aTol_crystalliteStress .or. &      ! below absolute tolerance .or.
         ( any(abs(dt*Lpguess) > relevantStrain) .and. &            ! worth checking? .and.
             maxval(abs(residuum/Lpguess), abs(dt*Lpguess) > relevantStrain) < rTol_crystalliteStress & ! below relative tolerance
         ) &
       ) &
      ) then
    exit LpLoop
  endif
  
  
  !* NaN occured at regular speed?  ->  return
  
  if (any(residuum/=residuum) .and. leapfrog == 1.0) then
#ifndef _OPENMP
    if (debug_verbosity > 4) then 
      write(6,'(a,i5,x,i2,x,i3,a,i3,a)') '<< CRYST >> integrateStress encountered NaN at el ip g ',e,i,g,&
                                    ' ; iteration ', NiterationStress,&
                                    ' >> returning..!'
    endif
#endif
    return

   
  !* something went wrong at accelerated speed?  ->  restore old residuum and Lp
 
  elseif (leapfrog > 1.0_pReal .and. &                             ! at fast pace .and.
           ( sum(residuum*residuum) > sum(residuum_old*residuum_old) .or. &  ! worse residuum .or.
             sum(residuum*residuum_old) < 0.0_pReal .or. &         ! residuum changed sign (overshoot) .or.
             any(residuum/=residuum) &                             ! NaN occured
           ) &
         ) then
#ifndef _OPENMP
    if (debug_verbosity > 5) then 
      write(6,'(a,i5,x,i2,x,i3,x,a,i3)') '<< CRYST >> integrateStress encountered high-speed crash at el ip g ',e,i,g,&
                                         '; iteration ', NiterationStress
    endif
#endif
    maxleap = 0.5_pReal * leapfrog                                 ! limit next acceleration
    leapfrog = 1.0_pReal                                           ! grinding halt
    jacoCounter = 0_pInt                                           ! reset counter for Jacobian update (we want to do an update next time!)
    Lpguess = Lpguess_old                                       
    residuum  = residuum_old
    if (debug_verbosity > 0) then
      !$OMP CRITICAL (distributionLeapfrogBreak)
        debug_LeapfrogBreakDistribution(NiterationStress,numerics_integrationMode) = &
          debug_LeapfrogBreakDistribution(NiterationStress,numerics_integrationMode) + 1
      !$OMP END CRITICAL (distributionLeapfrogBreak)
    endif
  
  
  !* residuum got better?  ->  calculate Jacobian for correction term  and remember current residuum and Lpguess
  
  else
    if (mod(jacoCounter, iJacoLpresiduum) == 0_pInt) then
      dTdLp = 0.0_pReal
      do h=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3
        dTdLp(3*(h-1)+j,3*(k-1)+l) = dTdLp(3*(h-1)+j,3*(k-1)+l) + C(h,j,l,m)*AB(k,m)+C(h,j,m,l)*BTA(m,k)
      enddo; enddo; enddo; enddo; enddo
      dTdLp = -0.5_pReal*dt*dTdLp
      dRdLp = math_identity2nd(9) - math_mul99x99(dLpdT_constitutive,dTdLp)
      invdRdLp = 0.0_pReal
      call math_invert(9,dRdLp,invdRdLp,dummy,error)               ! invert dR/dLp --> dLp/dR
      if (error) then
#ifndef _OPENMP
        if (debug_verbosity > 4) then
          write(6,'(a,i5,x,i2,x,i3,a,i3)') '<< CRYST >> integrateStress failed on dR/dLp inversion at el ip g ',e,i,g,&
                                      ' ; iteration ', NiterationStress
          if (debug_verbosity > 5 &
              .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
            write(6,*)
            write(6,'(a,/,9(12(x),9(e15.3,x)/))') '<< CRYST >> dRdLp',transpose(dRdLp)
            write(6,'(a,/,9(12(x),9(e15.3,x)/))') '<< CRYST >> dLpdT_constitutive',transpose(dLpdT_constitutive)
            write(6,'(a,/,3(12(x),3(e20.7,x)/))') '<< CRYST >> Lp_constitutive',math_transpose3x3(Lp_constitutive)
            write(6,'(a,/,3(12(x),3(e20.7,x)/))') '<< CRYST >> Lpguess',math_transpose3x3(Lpguess)
          endif
        endif
#endif
        return
      endif
    endif
    jacoCounter = jacoCounter + 1_pInt                             ! increase counter for jaco update
    residuum_old = residuum
    Lpguess_old = Lpguess 
     
    
    !* accelerate?
    
    if (NiterationStress > 1 .and. leapfrog+1.0_pReal <= maxleap) leapfrog = leapfrog + 1.0_pReal
    
  endif

  
  !* leapfrog to updated Lp
  
  do k=1,3; do l=1,3; do m=1,3; do n=1,3
    Lpguess(k,l) = Lpguess(k,l) - leapfrog * invdRdLp(3*(k-1)+l,3*(m-1)+n) * residuum(m,n)
  enddo; enddo; enddo; enddo

enddo LpLoop


!* calculate new plastic and elastic deformation gradient

invFp_new = math_mul33x33(invFp_current,B)
invFp_new = invFp_new/math_det3x3(invFp_new)**(1.0_pReal/3.0_pReal)  ! regularize by det
call math_invert3x3(invFp_new,Fp_new,det,error)
if (error) then
#ifndef _OPENMP
  if (debug_verbosity > 4) then
    write(6,'(a,i5,x,i2,x,i3,a,i3)') '<< CRYST >> integrateStress failed on invFp_new inversion at el ip g ',e,i,g, &
                                         ' ; iteration ', NiterationStress
    if (debug_verbosity > 5 .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger)) then
      write(6,*)
      write(6,'(a,/,3(12(x),3(f12.7,x)/))') '<< CRYST >> invFp_new',math_transpose3x3(invFp_new)
    endif
  endif
#endif
  return
endif
Fe_new = math_mul33x33(Fg_new,invFp_new)                             ! calc resulting Fe


!* add volumetric component to 2nd Piola-Kirchhoff stress and calculate 1st Piola-Kirchhoff stress

forall (n=1:3) Tstar_v(n) = Tstar_v(n) + p_hydro
crystallite_P(1:3,1:3,g,i,e) = math_mul33x33(Fe_new, math_mul33x33(math_Mandel6to33(Tstar_v), math_transpose3x3(invFp_new)))
 

!* store local values in global variables

crystallite_Lp(1:3,1:3,g,i,e) = Lpguess
crystallite_Tstar_v(1:6,g,i,e) = Tstar_v
crystallite_Fp(1:3,1:3,g,i,e) = Fp_new
crystallite_Fe(1:3,1:3,g,i,e) = Fe_new
crystallite_invFp(1:3,1:3,g,i,e) = invFp_new


!* set return flag to true

crystallite_integrateStress = .true.
#ifndef _OPENMP
if (debug_verbosity > 5 .and. ((e == debug_e .and. i == debug_i .and. g == debug_g) .or. .not. debug_selectiveDebugger) &
    .and. numerics_integrationMode == 1_pInt) then 
  write(6,'(a,/,3(12(x),3(f12.7,x)/))') '<< CRYST >> P / MPa',math_transpose3x3(crystallite_P(1:3,1:3,g,i,e))/1e6
  write(6,'(a,/,3(12(x),3(f12.7,x)/))') '<< CRYST >> Cauchy / MPa', &
                              math_mul33x33(crystallite_P(1:3,1:3,g,i,e), math_transpose3x3(Fg_new)) / 1e6 / math_det3x3(Fg_new)
  write(6,'(a,/,3(12(x),3(f12.7,x)/))') '<< CRYST >> Fe Lp Fe^-1', &
                      math_transpose3x3(math_mul33x33(Fe_new, math_mul33x33(crystallite_Lp(1:3,1:3,g,i,e), math_inv3x3(Fe_new))))                     ! transpose to get correct print out order
  write(6,'(a,/,3(12(x),3(f12.7,x)/))') '<< CRYST >> Fp',math_transpose3x3(crystallite_Fp(1:3,1:3,g,i,e))
endif
#endif

if (debug_verbosity > 0) then
  !$OMP CRITICAL (distributionStress)
   debug_StressLoopDistribution(NiterationStress,numerics_integrationMode) = &
     debug_StressLoopDistribution(NiterationStress,numerics_integrationMode) + 1
  !$OMP END CRITICAL (distributionStress)
endif

endfunction
 
 
 
!********************************************************************
! calculates orientations and disorientations (in case of single grain ips)
!******************************************************************** 
subroutine crystallite_orientations()
  
!*** variables and functions from other modules ***!
use prec, only:                       pInt, &
                                      pReal
use math, only:                       math_pDecomposition, &
                                      math_RtoQuaternion, &
                                      math_QuaternionDisorientation, &
                                      inDeg, &
                                      math_qConj
use FEsolving, only:                  FEsolving_execElem, & 
                                      FEsolving_execIP
use IO, only:                         IO_warning
use material, only:                   material_phase, &
                                      homogenization_Ngrains, &
                                      phase_constitution, &
                                      phase_localConstitution, &
                                      phase_constitutionInstance
use mesh, only:                       mesh_element, &
                                      mesh_ipNeighborhood, &
                                      FE_NipNeighbors
use debug, only:                      debug_verbosity, &
                                      debug_selectiveDebugger, &
                                      debug_e, debug_i, debug_g
use constitutive_nonlocal, only:      constitutive_nonlocal_structure, &
                                      constitutive_nonlocal_updateCompatibility

implicit none

!*** input variables ***!

!*** output variables ***!

!*** local variables ***!
integer(pInt)                   e, &                          ! element index
                                i, &                          ! integration point index
                                g, &                          ! grain index
                                n, &                          ! neighbor index 
                                neighboring_e, &              ! element index of my neighbor
                                neighboring_i, &              ! integration point index of my neighbor
                                myPhase, &                    ! phase
                                neighboringPhase, &
                                myInstance, &                 ! instance of constitution
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
      do g = 1,homogenization_Ngrains(mesh_element(3,e))
        
        call math_pDecomposition(crystallite_Fe(1:3,1:3,g,i,e), U, R, error)                              ! polar decomposition of Fe
        if (error) then
          call IO_warning(650, e, i, g)
          orientation = (/1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal/)                                    ! fake orientation
        else
          orientation = math_RtoQuaternion(transpose(R))
        endif
        crystallite_rotation(1:4,g,i,e) = math_QuaternionDisorientation(math_qConj(orientation), &
                                                                        math_qConj(crystallite_orientation0(1:4,g,i,e)), &
                                                                        0_pInt )                          ! we don't want symmetry here  
        crystallite_orientation(1:4,g,i,e) = orientation
      enddo
    enddo
  enddo
!$OMP END PARALLEL DO


! --- UPDATE SOME ADDITIONAL VARIABLES THAT ARE NEEDED FOR NONLOCAL MATERIAL ---
! --- we use crystallite_orientation from above, so need a seperate loop

!$OMP PARALLEL DO PRIVATE(myPhase,myInstance,myStructure,neighboring_e,neighboring_i, & 
!$OMP &                   neighboringPhase,neighboringInstance,neighboringStructure)
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      myPhase = material_phase(1,i,e)                                                                     ! get my phase
      if (.not. phase_localConstitution(myPhase)) then                                                    ! if nonlocal model
        myInstance = phase_constitutionInstance(myPhase)
        myStructure = constitutive_nonlocal_structure(myInstance)                                         ! get my crystal structure

        
        ! --- calculate disorientation between me and my neighbor ---
        
        do n = 1,FE_NipNeighbors(mesh_element(2,e))                                                       ! loop through my neighbors          
          neighboring_e = mesh_ipNeighborhood(1,n,i,e)
          neighboring_i = mesh_ipNeighborhood(2,n,i,e)
          if ((neighboring_e > 0) .and. (neighboring_i > 0)) then                                         ! if neighbor exists
            neighboringPhase = material_phase(1,neighboring_i,neighboring_e)                              ! get my neighbor's phase
            if (.not. phase_localConstitution(neighboringPhase)) then                                     ! neighbor got also nonlocal constitution
              neighboringInstance = phase_constitutionInstance(neighboringPhase)        
              neighboringStructure = constitutive_nonlocal_structure(neighboringInstance)                 ! get my neighbor's crystal structure               
              if (myStructure == neighboringStructure) then                                               ! if my neighbor has same crystal structure like me
                crystallite_disorientation(:,n,1,i,e) = &
                  math_QuaternionDisorientation( crystallite_orientation(1:4,1,i,e), &
                                                 crystallite_orientation(1:4,1,neighboring_i,neighboring_e), & 
                                                 crystallite_symmetryID(1,i,e))                           ! calculate disorientation            
              else                                                                                        ! for neighbor with different phase
                crystallite_disorientation(:,n,1,i,e) = (/0.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal/)    ! 180 degree rotation about 100 axis
              endif
            else                                                                                          ! for neighbor with local constitution
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

endsubroutine


 
!********************************************************************
! return results of particular grain
!********************************************************************
function crystallite_postResults(&
   dt,&             ! time increment
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )

 !*** variables and functions from other modules ***!
 use prec, only:                      pInt, &
                                      pReal
 use math, only:                      math_QuaternionToEuler, &
                                      math_QuaternionToAxisAngle, &
                                      math_mul33x33, &
                                      math_transpose3x3, &
                                      math_I3, &
                                      inDeg, &
                                      math_Mandel6to33
 use mesh, only:                      mesh_element
 use material, only:                  microstructure_crystallite, &
                                      crystallite_Noutput, &
                                      material_phase, &
                                      material_texture, &
                                      material_volume
 use constitutive, only:              constitutive_sizePostResults, &
                                      constitutive_postResults
 
 implicit none

 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 real(pReal), intent(in)::            dt                            ! time increment

 !*** output variables ***!
 real(pReal), dimension(1+crystallite_sizePostResults(microstructure_crystallite(mesh_element(4,e)))+ &
                        1+constitutive_sizePostResults(g,i,e)) :: crystallite_postResults
 
 !*** local variables ***!
 real(pReal), dimension(3,3) ::       Ee
 integer(pInt)                        o,c,crystID,mySize

 crystID = microstructure_crystallite(mesh_element(4,e))

 crystallite_postResults = 0.0_pReal
 c = 0_pInt
 crystallite_postResults(c+1) = crystallite_sizePostResults(crystID)            ! size of results from cryst
 c = c + 1_pInt
 
 do o = 1,crystallite_Noutput(crystID)
   select case(crystallite_output(o,crystID))
     case ('phase')
       crystallite_postResults(c+1) = material_phase(g,i,e)                     ! phaseID of grain
       c = c + 1_pInt
     case ('texture')
       crystallite_postResults(c+1) = material_texture(g,i,e)                   ! textureID of grain
       c = c + 1_pInt
     case ('volume')
       crystallite_postResults(c+1) = material_volume(g,i,e)                    ! grain volume (not fraction but absolute, right?)
       c = c + 1_pInt
     case ('orientation')
       crystallite_postResults(c+1:c+4) = crystallite_orientation(1:4,g,i,e)    ! grain orientation as quaternion
       c = c + 4_pInt
     case ('eulerangles')
       crystallite_postResults(c+1:c+3) = inDeg * math_QuaternionToEuler(crystallite_orientation(1:4,g,i,e)) ! grain orientation as Euler angles in degree
       c = c + 3_pInt
     case ('grainrotation')
       crystallite_postResults(c+1:c+4) = math_QuaternionToAxisAngle(crystallite_rotation(1:4,g,i,e)) ! grain rotation away from initial orientation as axis-angle 
       crystallite_postResults(c+4) = inDeg * crystallite_postResults(c+4)      ! angle in degree
       c = c + 4_pInt
     case ('defgrad','f')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_transpose3x3(crystallite_partionedF(1:3,1:3,g,i,e)),(/mySize/))
       c = c + mySize
     case ('fe')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_transpose3x3(crystallite_Fe(1:3,1:3,g,i,e)),(/mySize/))
       c = c + mySize
     case ('ee')
       Ee = 0.5_pReal * (math_mul33x33(math_transpose3x3(crystallite_Fe(1:3,1:3,g,i,e)), crystallite_Fe(1:3,1:3,g,i,e)) - math_I3)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(Ee,(/mySize/))
       c = c + mySize
     case ('fp')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_transpose3x3(crystallite_Fp(1:3,1:3,g,i,e)),(/mySize/))
       c = c + mySize
     case ('lp')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_transpose3x3(crystallite_Lp(1:3,1:3,g,i,e)),(/mySize/))
       c = c + mySize
     case ('p','firstpiola','1stpiola')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_transpose3x3(crystallite_P(1:3,1:3,g,i,e)),(/mySize/))
       c = c + mySize
     case ('s','tstar','secondpiola','2ndpiola')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_Mandel6to33(crystallite_Tstar_v(1:6,g,i,e)),(/mySize/))
       c = c + mySize
   end select
 enddo
  
 crystallite_postResults(c+1) = constitutive_sizePostResults(g,i,e)             ! size of constitutive results
 c = c + 1_pInt
 crystallite_postResults(c+1:c+constitutive_sizePostResults(g,i,e)) = constitutive_postResults(crystallite_Tstar_v(1:6,g,i,e), &
                                                                                               crystallite_Fe(1:3,1:3,g,i,e), &
                                                                                               crystallite_Temperature(g,i,e), &
                                                                                               dt, g, i, e)
 c = c + constitutive_sizePostResults(g,i,e)
 
 return 
 
endfunction


END MODULE
!##############################################################

