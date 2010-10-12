!* $Id$
!***************************************
!*      Module: CRYSTALLITE            *
!***************************************
!* contains:                           *
!* - _init                             *
!* - materialpoint_stressAndItsTangent *
!* - _partitionDeformation             *
!* - _updateState                      *
!* - _stressAndItsTangent       *
!* - _postResults                      *
!***************************************

MODULE crystallite

use prec, only: pReal, pInt
implicit none
!
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
    crystallite_converged, &             ! convergence flag
    crystallite_stateConverged, &        ! flag indicating convergence of state
    crystallite_temperatureConverged     ! flag indicating convergence of temperature

CONTAINS

!********************************************************************
! allocate and initialize per grain variables
!********************************************************************
subroutine crystallite_init(Temperature)
  
 !*** variables and functions from other modules ***!
 use prec, only:             pInt, &
                             pReal
 use debug, only:            debug_info, &
                             debug_reset
 use numerics, only:         integrator, &
                             integratorStiffness, &
                             subStepSizeCryst, &
                             stepIncreaseCryst
 use math, only:             math_I3, &
                             math_EulerToR
 use FEsolving, only:        FEsolving_execElem, &
                             FEsolving_execIP
 use mesh, only:             mesh_element, &
                             mesh_NcpElems, &
                             mesh_maxNips, &
                             mesh_maxNipNeighbors
 use IO
 use material
 use lattice, only:          lattice_symmetryTypes
 use constitutive_phenopowerlaw, only: constitutive_phenopowerlaw_label, &
                                      constitutive_phenopowerlaw_structure
 use constitutive_titanmod, only: constitutive_titanmod_label, &
                                      constitutive_titanmod_structure
 use constitutive_dislotwin, only:     constitutive_dislotwin_label, &
                                      constitutive_dislotwin_structure
 use constitutive_nonlocal, only:      constitutive_nonlocal_label, &
                                      constitutive_nonlocal_structure
 
 implicit none
 integer(pInt), parameter :: file = 200
 
 !*** input variables ***!
 real(pReal) Temperature
 
 !*** output variables ***!

 !*** local variables ***!
 integer(pInt), parameter :: maxNchunks = 2
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt)               g, &                          ! grain number
                             i, &                          ! integration point number
                             e, &                          ! element number
                             gMax, &                       ! maximum number of grains
                             iMax, &                       ! maximum number of integration points
                             eMax, &                       ! maximum number of elements
                             nMax, &                       ! maximum number of ip neighbors
                             myNgrains, &                  ! number of grains in current IP
                             myCrystallite                 ! crystallite of current elem
 integer(pInt) section, j,p, output, mySize
 character(len=64) tag
 character(len=1024) line
 integer(pInt)               myStructure, &                ! lattice structure 
                             myPhase
 
 write(6,*)
 write(6,*) '<<<+-  crystallite init  -+>>>'
 write(6,*) '$Id$'
 write(6,*)
 
 
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
 allocate(crystallite_stateConverged(gMax,iMax,eMax));               crystallite_stateConverged = .false.
 allocate(crystallite_temperatureConverged(gMax,iMax,eMax));   crystallite_temperatureConverged = .false.
 
 allocate(crystallite_output(maxval(crystallite_Noutput), &
                             material_Ncrystallite)) ;                       crystallite_output = ''
 allocate(crystallite_sizePostResults(material_Ncrystallite)) ;     crystallite_sizePostResults = 0_pInt
 allocate(crystallite_sizePostResult(maxval(crystallite_Noutput), &
                             material_Ncrystallite)) ;               crystallite_sizePostResult = 0_pInt
 
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
       case('phase')
         mySize = 1
       case('volume')
         mySize = 1
       case('orientation')   ! orientation as quaternion
         mySize = 4
       case('eulerangles')   ! Bunge Euler angles
         mySize = 3
       case('grainrotation') ! Deviation from initial grain orientation in axis-angle form (angle in degrees)
         mySize = 4
       case('defgrad','f','fe','fp','ee','p','firstpiola','1stpiola','s','tstar','secondpiola','2ndpiola')
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


 !$OMP PARALLEL DO
 do e = FEsolving_execElem(1),FEsolving_execElem(2)                       ! iterate over all cp elements
   myNgrains = homogenization_Ngrains(mesh_element(3,e))                  ! look up homogenization-->grainCount
   do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                     ! iterate over IPs of this element
     do g = 1,myNgrains
       crystallite_partionedTemperature0(g,i,e) = Temperature                       ! isothermal assumption
       crystallite_Fp0(:,:,g,i,e) = math_EulerToR(material_EulerAngles(:,g,i,e))    ! plastic def gradient reflects init orientation
       crystallite_Fe(:,:,g,i,e)  = transpose(crystallite_Fp0(:,:,g,i,e))
       crystallite_F0(:,:,g,i,e)  = math_I3
       crystallite_partionedFp0(:,:,g,i,e) = crystallite_Fp0(:,:,g,i,e)
       crystallite_partionedF0(:,:,g,i,e)  = crystallite_F0(:,:,g,i,e)
       crystallite_partionedF(:,:,g,i,e)   = crystallite_F0(:,:,g,i,e)
       crystallite_requested(g,i,e)        = .true.
       crystallite_localConstitution(g,i,e) = phase_localConstitution(material_phase(g,i,e))
     enddo
   enddo
 enddo
 !$OMPEND PARALLEL DO

! Initialize crystallite_symmetryID(g,i,e)
!$OMP PARALLEL DO 
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do g = 1,homogenization_Ngrains(mesh_element(3,e))
         myPhase = material_phase(g,i,e)
         select case (phase_constitution(myPhase))
                case (constitutive_phenopowerlaw_label)
                  myStructure = constitutive_phenopowerlaw_structure(phase_constitutionInstance(myPhase))
                case (constitutive_titanmod_label)
                  myStructure = constitutive_titanmod_structure(phase_constitutionInstance(myPhase))
                case (constitutive_dislotwin_label)
                  myStructure = constitutive_dislotwin_structure(phase_constitutionInstance(myPhase))
                case (constitutive_nonlocal_label)
                  myStructure = constitutive_nonlocal_structure(phase_constitutionInstance(myPhase))
                case default
                  myStructure = -1_pInt ! does this happen for j2 material?
              end select
         if (myStructure>0_pInt) then   
           crystallite_symmetryID(g,i,e)=lattice_symmetryTypes(myStructure) ! structure = 1(fcc) or 2(bcc) => 1; 3(hex)=>2  
         endif
      enddo
    enddo
  enddo
!$OMPEND PARALLEL DO         
 
  call crystallite_orientations()
  crystallite_orientation0 = crystallite_orientation  ! Store initial orientations for calculation of grain rotations

  call crystallite_stressAndItsTangent(.true.)                 ! request elastic answers
  crystallite_fallbackdPdF = crystallite_dPdF                  ! use initial elastic stiffness as fallback
 
  !    *** Output to MARC output file ***
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
    write(6,'(a35,x,7(i5,x))') 'crystallite_stateConverged:        ', shape(crystallite_stateConverged)
    write(6,'(a35,x,7(i5,x))') 'crystallite_temperatureConverged:  ', shape(crystallite_temperatureConverged)
    write(6,'(a35,x,7(i5,x))') 'crystallite_sizePostResults:       ', shape(crystallite_sizePostResults)
    write(6,'(a35,x,7(i5,x))') 'crystallite_sizePostResult:        ', shape(crystallite_sizePostResult)
    write(6,*)
    write(6,*) 'Number of nonlocal grains: ',count(.not. crystallite_localConstitution)
    call flush(6)
  !$OMPEND CRITICAL (write2out)

  call debug_info()
  call debug_reset()

  return

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
                                                        integratorStiffness, &
                                                        integrator
  use debug, only:                                      debugger, &
                                                        selectiveDebugger, &
                                                        verboseDebugger, &
                                                        debug_e, &
                                                        debug_i, &
                                                        debug_g, &
                                                        debug_CrystalliteLoopDistribution
  use IO, only:                                         IO_warning
  use math, only:                                       math_inv3x3, &
                                                        math_mul33x33, &
                                                        math_mul66x6, &
                                                        math_Mandel6to33, &
                                                        math_Mandel33to6, &
                                                        math_I3
  use FEsolving, only:                                  FEsolving_execElem, & 
                                                        FEsolving_execIP, &
                                                        theInc, &
                                                        cycleCounter
  use mesh, only:                                       mesh_element, &
                                                        mesh_NcpElems, &
                                                        mesh_maxNips
  use material, only:                                   homogenization_Ngrains, &
                                                        homogenization_maxNgrains
  use constitutive, only:                               constitutive_maxSizeState, &
                                                        constitutive_maxSizeDotState, &
                                                        constitutive_sizeState, &
                                                        constitutive_sizeDotState, &
                                                        constitutive_state, &
                                                        constitutive_state_backup, &
                                                        constitutive_subState0, &
                                                        constitutive_partionedState0, &
                                                        constitutive_homogenizedC, &
                                                        constitutive_dotState, &
                                                        constitutive_dotState_backup, &
                                                        constitutive_collectDotState, &
                                                        constitutive_dotTemperature, &
                                                        constitutive_microstructure

  implicit none

  !*** input variables ***!
  logical, intent(in) ::                                updateJaco                    ! flag indicating wehther we want to update the Jacobian (stiffness) or not

  !*** output variables ***!

  !*** local variables ***!
  real(pReal)                                           myTemperature, &              ! local copy of the temperature
                                                        myPert, &                     ! perturbation with correct sign
                                                        formerSubStep
  real(pReal), dimension(3,3) ::                        invFp, &                      ! inverse of the plastic deformation gradient
                                                        Fe_guess, &                   ! guess for elastic deformation gradient
                                                        Tstar                         ! 2nd Piola-Kirchhoff stress tensor
  real(pReal), dimension(9,9) ::                        dPdF99
  real(pReal), dimension(3,3,3,3,2) ::                  dPdF_perturbation
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
  logical                                               forceLocalStiffnessCalculation ! flag indicating that stiffness calculation is always done locally
  forceLocalStiffnessCalculation = .true. 

  ! --+>> INITIALIZE TO STARTING CONDITION <<+--
  
  crystallite_subStep = 0.0_pReal

  !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          if (crystallite_requested(g,i,e)) then                                                          ! initialize restoration point of ...
            crystallite_subTemperature0(g,i,e) = crystallite_partionedTemperature0(g,i,e)                 ! ...temperature
            constitutive_subState0(g,i,e)%p   = constitutive_partionedState0(g,i,e)%p                     ! ...microstructure
            crystallite_subFp0(:,:,g,i,e)     = crystallite_partionedFp0(:,:,g,i,e)                       ! ...plastic def grad
            crystallite_subLp0(:,:,g,i,e)     = crystallite_partionedLp0(:,:,g,i,e)                       ! ...plastic velocity grad
            crystallite_dPdF0(:,:,:,:,g,i,e)  = crystallite_partioneddPdF0(:,:,:,:,g,i,e)                 ! ...stiffness
            crystallite_subF0(:,:,g,i,e)      = crystallite_partionedF0(:,:,g,i,e)                        ! ...def grad
            crystallite_subTstar0_v(:,g,i,e)  = crystallite_partionedTstar0_v(:,g,i,e)                    !...2nd PK stress

            crystallite_subFrac(g,i,e) = 0.0_pReal
            crystallite_subStep(g,i,e) = 1.0_pReal/subStepSizeCryst
            crystallite_todo(g,i,e) = .true.
            crystallite_converged(g,i,e) = .false.                                                        ! pretend failed step of twice the required size
          endif
        enddo
      enddo
    enddo
  !$OMPEND PARALLEL DO
 

  ! --+>> CRYSTALLITE CUTBACK LOOP <<+--

  NiterationCrystallite = 0_pInt
  do while (any(crystallite_subStep(:,:,FEsolving_execELem(1):FEsolving_execElem(2)) > subStepMinCryst))  ! cutback loop for crystallites

    !$OMP PARALLEL DO
      do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                  ! iterate over elements to be processed
        myNgrains = homogenization_Ngrains(mesh_element(3,e))
        do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                ! iterate over IPs of this element to be processed
          do g = 1,myNgrains
            selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g) 
            
            ! --- wind forward ---
            
            if (crystallite_converged(g,i,e)) then
              if (debugger .and. selectiveDebugger) then
                !$OMP CRITICAL (write2out)
                  write(6,'(a21,f10.8,a32,f10.8,a35)') 'winding forward from ', &
                    crystallite_subFrac(g,i,e),' to current crystallite_subfrac ', &
                    crystallite_subFrac(g,i,e)+crystallite_subStep(g,i,e),' in crystallite_stressAndItsTangent'
                  write(6,*)
                !$OMPEND CRITICAL (write2out)
              endif
              crystallite_subFrac(g,i,e) = crystallite_subFrac(g,i,e) + crystallite_subStep(g,i,e)
              formerSubStep = crystallite_subStep(g,i,e)
              crystallite_subStep(g,i,e) = min( 1.0_pReal - crystallite_subFrac(g,i,e), &
                                                stepIncreaseCryst * crystallite_subStep(g,i,e) )
              if (crystallite_subStep(g,i,e) > subStepMinCryst) then
                crystallite_subTemperature0(g,i,e) = crystallite_Temperature(g,i,e)                       ! wind forward...
                crystallite_subF0(:,:,g,i,e) = crystallite_subF(:,:,g,i,e)                                ! ...def grad
                crystallite_subFp0(:,:,g,i,e) = crystallite_Fp(:,:,g,i,e)                                 ! ...plastic def grad
                crystallite_subLp0(:,:,g,i,e) = crystallite_Lp(:,:,g,i,e)                                 ! ...plastic velocity gradient
                constitutive_subState0(g,i,e)%p = constitutive_state(g,i,e)%p                             ! ...microstructure
                crystallite_subTstar0_v(:,g,i,e) = crystallite_Tstar_v(:,g,i,e)                           ! ...2nd PK stress
              elseif (formerSubStep > subStepMinCryst) then                                               ! this crystallite just converged
               !$OMP CRITICAL (distributionCrystallite)
                 debug_CrystalliteLoopDistribution(min(nCryst+1,NiterationCrystallite)) = &
                   debug_CrystalliteLoopDistribution(min(nCryst+1,NiterationCrystallite)) + 1
               !$OMPEND CRITICAL (distributionCrystallite)
              endif
            
            ! --- cutback ---
            
            else
              crystallite_subStep(g,i,e) = subStepSizeCryst * crystallite_subStep(g,i,e)                  ! cut step in half and restore...
              crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e)                         ! ...temperature
              crystallite_Fp(:,:,g,i,e) = crystallite_subFp0(:,:,g,i,e)                                   ! ...plastic def grad
              crystallite_invFp(:,:,g,i,e) = math_inv3x3(crystallite_Fp(:,:,g,i,e))
              crystallite_Lp(:,:,g,i,e) = crystallite_subLp0(:,:,g,i,e)                                   ! ...plastic velocity grad
              constitutive_state(g,i,e)%p = constitutive_subState0(g,i,e)%p                               ! ...microstructure
              crystallite_Tstar_v(:,g,i,e) = crystallite_subTstar0_v(:,g,i,e)                             ! ...2nd PK stress
                                                                                                          ! can�t restore dotState here, since not yet calculated in first cutback after initialization
              if (debugger .and. selectiveDebugger) then
                !$OMP CRITICAL (write2out)
                  write(6,'(a78,f10.8)') 'cutback step in crystallite_stressAndItsTangent with new crystallite_subStep: ',&
                                         crystallite_subStep(g,i,e)
                  write(6,*)
                !$OMPEND CRITICAL (write2out)
              endif
            endif

            ! --- prepare for integration ---
            
            crystallite_todo(g,i,e) = crystallite_subStep(g,i,e) > subStepMinCryst                     ! still on track or already done (beyond repair)
            if (crystallite_todo(g,i,e)) then
              crystallite_subF(:,:,g,i,e)  = crystallite_subF0(:,:,g,i,e) + &
                                             crystallite_subStep(g,i,e) * &
                                             (crystallite_partionedF(:,:,g,i,e) - crystallite_partionedF0(:,:,g,i,e))
              crystallite_Fe(:,:,g,i,e)    = math_mul33x33(crystallite_subF(:,:,g,i,e),crystallite_invFp(:,:,g,i,e))
              crystallite_subdt(g,i,e)     = crystallite_subStep(g,i,e) * crystallite_dt(g,i,e)
              crystallite_converged(g,i,e) = .false.                                                      ! start out non-converged
            endif
            
          enddo
        enddo
      enddo
    !$OMPEND PARALLEL DO

    ! --- integrate ---
    
    if (any(crystallite_todo)) then
      select case(integrator)
        case (1)
          call crystallite_integrateStateFPI(1)
        case (2)
          call crystallite_integrateStateEuler(1)
        case (3)
          call crystallite_integrateStateAdaptiveEuler(1)
        case (4)
          call crystallite_integrateStateRK4(1)
        case(5)
          call crystallite_integrateStateRKCK45(1)
      endselect
    endif
    
    NiterationCrystallite = NiterationCrystallite + 1
        
  enddo                                                                                                   ! cutback loop
  
  
  ! --+>> CHECK FOR NON-CONVERGED CRYSTALLITES <<+--

  !$OMP PARALLEL DO
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          if (.not. crystallite_converged(g,i,e)) then                                                    ! respond fully elastically (might be not required due to becoming terminally ill anyway)
!            call IO_warning(600,e,i,g)
            invFp = math_inv3x3(crystallite_partionedFp0(:,:,g,i,e))
            Fe_guess = math_mul33x33(crystallite_partionedF(:,:,g,i,e),invFp)
            Tstar = math_Mandel6to33( math_mul66x6( 0.5_pReal*constitutive_homogenizedC(g,i,e), &
                                                    math_Mandel33to6( math_mul33x33(transpose(Fe_guess),Fe_guess) - math_I3 ) ) )
            crystallite_P(:,:,g,i,e) = math_mul33x33(Fe_guess,math_mul33x33(Tstar,transpose(invFp)))
          endif
        enddo
      enddo
    enddo
  !$OMPEND PARALLEL DO


  ! --+>> STIFFNESS CALCULATION <<+--
  
  if(updateJaco) then                                                                                     ! Jacobian required
    
    ! --- BACKUP ---
    
    do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                    ! iterate over elements to be processed
      myNgrains = homogenization_Ngrains(mesh_element(3,e))
      do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                                  ! iterate over IPs of this element to be processed
        do g = 1,myNgrains
          mySizeState = constitutive_sizeState(g,i,e)                                                     ! number of state variables for this grain
          mySizeDotState = constitutive_sizeDotState(g,i,e)                                               ! number of dotStates for this grain
          constitutive_state_backup(g,i,e)%p(1:mySizeState) = &
                 constitutive_state(g,i,e)%p(1:mySizeState)                                               ! remember unperturbed, converged state, ...
          constitutive_dotState_backup(g,i,e)%p(1:mySizeDotState) = &
                 constitutive_dotState(g,i,e)%p(1:mySizeDotState)                                         ! ... dotStates, ...
    enddo; enddo; enddo
    Temperature_backup = crystallite_Temperature                                                          ! ... Temperature, ...
    F_backup = crystallite_subF                                                                           ! ... and kinematics
    Fp_backup = crystallite_Fp
    InvFp_backup = crystallite_invFp
    Fe_backup = crystallite_Fe
    Lp_backup = crystallite_Lp
    Tstar_v_backup = crystallite_Tstar_v
    P_backup = crystallite_P
    convergenceFlag_backup = crystallite_converged

    
    ! --- LOCAL STIFFNESS CALCULATION ---
    
    if (all(crystallite_localConstitution) .or. theInc < 1 .or. forceLocalStiffnessCalculation) then      ! all grains have local constitution, so local convergence of perturbed grain is sufficient
    
      !$OMP PARALLEL DO
        do e = FEsolving_execElem(1),FEsolving_execElem(2)                                                ! iterate over elements to be processed
          myNgrains = homogenization_Ngrains(mesh_element(3,e))
          do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)                                              ! iterate over IPs of this element to be processed
            do g = 1,myNgrains
              selectiveDebugger = .false. ! (e == debug_e .and. i == debug_i .and. g == debug_g)
              if (crystallite_requested(g,i,e)) then                                                      ! first check whether is requested at all!
                if (crystallite_converged(g,i,e)) then                                                    ! grain converged in above iteration
                  
                  if (verboseDebugger .and. selectiveDebugger) then
                    !$OMP CRITICAL (write2out)
                    write (6,*) '#############'
                    write (6,*) 'central solution of cryst_StressAndTangent'
                    write (6,*) '#############'
                    write (6,'(a8,3(x,i4),/,3(3(f12.4,x)/))') '    P of', g, i, e, P_backup(1:3,:,g,i,e)/1e6
                    write (6,'(a8,3(x,i4),/,3(3(f14.9,x)/))') '   Fp of', g, i, e, Fp_backup(1:3,:,g,i,e)
                    write (6,'(a8,3(x,i4),/,3(3(f14.9,x)/))') '   Lp of', g, i, e, Lp_backup(1:3,:,g,i,e)
                    !$OMPEND CRITICAL (write2out)
                  endif
                  
                  do perturbation = 1,2                                                                     ! forward and backward perturbation
                    if (iand(pert_method,perturbation) > 0) then                                            ! mask for desired direction
                      dPdF_perturbation(:,:,:,:,perturbation) = crystallite_dPdF0(:,:,:,:,g,i,e)            ! initialize stiffness with known good values from last increment
                      myPert = -pert_Fg * (-1.0_pReal)**perturbation                                        ! set perturbation step
                      do k = 1,3; do l = 1,3                                                                ! ...alter individual components
                        crystallite_subF(k,l,g,i,e) = crystallite_subF(k,l,g,i,e) + myPert                  ! perturb either forward or backward

                        if (verboseDebugger .and. selectiveDebugger) then
                          !$OMP CRITICAL (write2out)
                          write (6,'(a,x,i1,x,i1,x,a)') '[[[[[[[ Stiffness perturbation',k,l,']]]]]]]'
                          write (6,'(a8,3(x,i4),/,3(3(f14.9,x)/))') 'pertF of', g, i, e, crystallite_subF(1:3,:,g,i,e)
                          !$OMPEND CRITICAL (write2out)
                        endif
                        
                        ! --- local integration and stiffness calculation ---
                        
                        crystallite_converged(g,i,e) = .false.                                                                    ! start out non-converged
                        crystallite_todo(g,i,e) = .true.
                        select case(integratorStiffness)
                          case (1)
                            call crystallite_integrateStateFPI(2,g,i,e)
                          case (2)
                            call crystallite_integrateStateEuler(2,g,i,e)
                          case (3)
                            call crystallite_integrateStateAdaptiveEuler(2,g,i,e)
                          case (4)
                            call crystallite_integrateStateRK4(2,g,i,e)
                          case(5)
                            call crystallite_integrateStateRKCK45(2,g,i,e)
                        endselect
                        if (crystallite_converged(g,i,e)) &                                                                       ! converged state warrants stiffness update
                          dPdF_perturbation(:,:,k,l,perturbation) = (crystallite_P(:,:,g,i,e) - P_backup(:,:,g,i,e))/myPert       ! tangent dP_ij/dFg_kl
  
                        ! --- restore ---
  
                        mySizeState = constitutive_sizeState(g,i,e)                                                               ! number of state variables for this grain
                        mySizeDotState = constitutive_sizeDotState(g,i,e)                                                         ! number of dotStates for this grain
                        constitutive_state(g,i,e)%p(1:mySizeState) = constitutive_state_backup(g,i,e)%p(1:mySizeState)
                        constitutive_dotState(g,i,e)%p(1:mySizeDotState) = constitutive_dotState_backup(g,i,e)%p(1:mySizeDotState)
                        crystallite_Temperature(g,i,e) = Temperature_backup(g,i,e)
                        crystallite_subF(:,:,g,i,e) = F_backup(:,:,g,i,e)
                        crystallite_Fp(:,:,g,i,e) = Fp_backup(:,:,g,i,e) 
                        crystallite_invFp(:,:,g,i,e) = InvFp_backup(:,:,g,i,e)
                        crystallite_Fe(:,:,g,i,e) = Fe_backup(:,:,g,i,e)
                        crystallite_Lp(:,:,g,i,e) = Lp_backup(:,:,g,i,e)
                        crystallite_Tstar_v(:,g,i,e) = Tstar_v_backup(:,g,i,e)
                        crystallite_P(:,:,g,i,e) = P_backup(:,:,g,i,e)
                        crystallite_converged(g,i,e) = convergenceFlag_backup(g,i,e)
                        
                      enddo; enddo
                    endif
                  enddo                                                                                 ! perturbation direction
                  select case(pert_method)
                    case (1)
                      crystallite_dPdF(:,:,:,:,g,i,e) = dPdF_perturbation(:,:,:,:,1)
                    case (2)
                      crystallite_dPdF(:,:,:,:,g,i,e) = dPdF_perturbation(:,:,:,:,2)
                    case (3)
                      crystallite_dPdF(:,:,:,:,g,i,e) = 0.5_pReal*(dPdF_perturbation(:,:,:,:,1)+dPdF_perturbation(:,:,:,:,2))
                  end select
                else                                                                                    ! grain did not converge
                  crystallite_dPdF(:,:,:,:,g,i,e) = crystallite_fallbackdPdF(:,:,:,:,g,i,e)             ! use (elastic) fallback
                endif               ! grain convergence
              endif                 ! grain request
            enddo                   ! grain   loop
          enddo                     ! ip      loop
        enddo                       ! element loop
      !$OMPEND PARALLEL DO
      

    ! --- NON-LOCAL STIFFNESS CALCULATION ---

    elseif (any(.not. crystallite_localConstitution)) then                                                                          ! if any nonlocal grain present, we have to do a full loop over all grains after each perturbance
      
      crystallite_dPdF = crystallite_dPdF0                                                                                          ! initialize stiffness with known good values from last inc

      do k = 1,3
        do l = 1,3
          crystallite_subF(k,l,:,:,:) = crystallite_subF(k,l,:,:,:) + pert_Fg                                                       ! perturb single component
                
          ! --- integration ---
          
          crystallite_converged = .false.                                                                                           ! start out non-converged
          crystallite_todo = .true.
          select case(integratorStiffness)
            case (1)
              call crystallite_integrateStateFPI(2)
            case (2)
              call crystallite_integrateStateEuler(2)
            case (3)
              call crystallite_integrateStateAdaptiveEuler(2)
            case (4)
              call crystallite_integrateStateRK4(2)
            case(5)
              call crystallite_integrateStateRKCK45(2)
          endselect
          
          ! --- stiffness calculation ---
          
          do e = FEsolving_execElem(1),FEsolving_execElem(2)
            myNgrains = homogenization_Ngrains(mesh_element(3,e))
            do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
              do g = 1,myNgrains
                if (crystallite_converged(g,i,e)) then                                                                              ! if stiffness calculation converged...
                  crystallite_dPdF(:,:,k,l,g,i,e) = (crystallite_P(:,:,g,i,e) - P_backup(:,:,g,i,e))/pert_Fg                        ! ... use tangent dP_ij/dFg_kl
                elseif (.not. convergenceFlag_backup(g,i,e)) then                                                                   ! if crystallite didnt converge before...
                  crystallite_dPdF(:,:,:,:,g,i,e) = crystallite_fallbackdPdF(:,:,:,:,g,i,e)                                         ! ... use (elastic) fallback
                endif
          enddo; enddo; enddo
          
          ! --- restore ---
          
          do e = FEsolving_execElem(1),FEsolving_execElem(2)
            myNgrains = homogenization_Ngrains(mesh_element(3,e))
            do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
              do g = 1,myNgrains
                mySizeState = constitutive_sizeState(g,i,e)
                mySizeDotState = constitutive_sizeDotState(g,i,e)
                constitutive_state(g,i,e)%p(1:mySizeState) = constitutive_state_backup(g,i,e)%p(1:mySizeState)
                constitutive_dotState(g,i,e)%p(1:mySizeDotState) = constitutive_dotState_backup(g,i,e)%p(1:mySizeDotState)
          enddo; enddo; enddo
          crystallite_Temperature = Temperature_backup
          crystallite_subF = F_backup
          crystallite_Fp = Fp_backup 
          crystallite_invFp = InvFp_backup
          crystallite_Fe = Fe_backup
          crystallite_Lp = Lp_backup
          crystallite_Tstar_v = Tstar_v_backup
          crystallite_P = P_backup
          crystallite_converged = convergenceFlag_backup
      
      enddo;enddo                 ! k,l loop
      
    endif

  endif                           ! jacobian calculation
 
endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 4h order explicit Runge Kutta method 
!********************************************************************
subroutine crystallite_integrateStateRK4(mode,gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debugger, &
                        selectiveDebugger, &
                        verboseDebugger, &
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
                        constitutive_dotTemperature, &
                        constitutive_microstructure

implicit none

real(pReal), dimension(4), parameter ::       timeStepFraction = (/0.5_pReal, 0.5_pReal, 1.0_pReal, 1.0_pReal/) ! weight of slope used for Runge Kutta integration
real(pReal), dimension(3), parameter ::       weight = (/2.0_pReal, 2.0_pReal, 1.0_pReal/)                      ! factor giving the fraction of the original timestep used for Runge Kutta Integration

!*** input variables ***!
integer(pInt), intent(in) ::                  mode                      ! mode of calculation; 1: central solution, 2: stiffness (by perturbation)
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
  iIter(:,ee) = ii
  gIter(:,ee) = gg
  singleRun = .true.
else
  eIter = FEsolving_execElem(1:2)
  do e = eIter(1),eIter(2)
    iIter(:,e) = FEsolving_execIP(1:2,e)
    gIter(:,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- UPDATE DEPENDENT STATES ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), crystallite_Fe, &
                                       crystallite_Fp, g, i, e)                                           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- FIRST RUNGE KUTTA STEP ---

RK4dotTemperature = 0.0_pReal
!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
      call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                        crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
      if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      else                                                                                                ! everything is fine
        constitutive_RK4dotState(g,i,e)%p = constitutive_dotState(g,i,e)%p                                ! initial contribution to RK slope
        RK4dotTemperature(g,i,e) = crystallite_dotTemperature(g,i,e)
      endif
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- SECOND TO FOURTH RUNGE KUTTA STEP PLUS FINAL INTEGRATION ---

do n = 1,4


  ! --- state update ---

  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        mySizeDotState = constitutive_sizeDotState(g,i,e)
        constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e) * timeStepFraction(n)
        crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e) * timeStepFraction(n)
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO      

  
  ! --- update dependent states ---

  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                         crystallite_Fe, crystallite_Fp, g, i, e)                         ! update dependent state variables to be consistent with basic states
      endif
   enddo; enddo; enddo
  !$OMPEND PARALLEL DO


  ! --- stress integration ---

  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
        if (crystallite_integrateStress(mode,g,i,e,timeStepFraction(n))) then                             ! fraction of original times step
          if (n == 4) then                                                                                ! final integration step
            if (mode==1 .and. verboseDebugger .and. e == debug_e .and. i == debug_i .and. g == debug_g) then
              mySizeDotState = constitutive_sizeDotState(g,i,e)
              !$OMP CRITICAL (write2out)
                write(6,*) '::: updateState',g,i,e
                write(6,*)
                write(6,'(a,/,12(e14.8,x))') 'updateState: dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
                write(6,*)
                write(6,'(a,/,12(e14.8,x))') 'updateState: new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
                write(6,*)
              !$OMPEND CRITICAL (write2out)
            endif
            crystallite_converged(g,i,e) = .true.                                                         ! ... converged per definition
            crystallite_todo(g,i,e) = .false.                                                             ! ... integration done
            !$OMP CRITICAL (distributionState)
              debug_StateLoopDistribution(n,mode) = debug_StateLoopDistribution(n,mode) + 1
            !$OMPEND CRITICAL (distributionState)
          endif
        else                                                                                              ! broken stress integration
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMPEND CRITICAL
          else                                                                                            ! if broken local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
          endif
        endif
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO      

  
  ! --- dot state and RK dot state---

  if (n < 4) then
    !$OMP PARALLEL DO
      do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                      ! iterate over elements, ips and grains
        if (crystallite_todo(g,i,e)) then
          selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
          call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), &
                                            crystallite_Fe, crystallite_Fp, crystallite_Temperature(g,i,e), & 
                                            timeStepFraction(n)*crystallite_subdt(g,i,e), g,i,e)          ! fraction of original timestep
          crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                          crystallite_Temperature(g,i,e),g,i,e)
          if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                      ! NaN occured in dotState
               .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then           ! NaN occured in dotTemperature
            if (.not. crystallite_localConstitution(g,i,e)) then                                          ! if broken non-local...
              !$OMP CRITICAL
                crystallite_todo = crystallite_todo .and. crystallite_localConstitution                   ! ...all non-locals skipped
              !$OMPEND CRITICAL
            else                                                                                          ! if broken local...
              crystallite_todo(g,i,e) = .false.                                                           ! ... skip this one next time
            endif
          else                                                                                            ! everything is fine
            constitutive_RK4dotState(g,i,e)%p = constitutive_RK4dotState(g,i,e)%p + weight(n)*constitutive_dotState(g,i,e)%p
            RK4dotTemperature(g,i,e) = RK4dotTemperature(g,i,e) + weight(n)*crystallite_dotTemperature(g,i,e)
            if (n == 3) then
              constitutive_dotState(g,i,e)%p = constitutive_RK4dotState(g,i,e)%p / 6.0_pReal              ! use weighted RKdotState for final integration
            endif
          endif
        endif
      enddo; enddo; enddo
    !$OMPEND PARALLEL DO
  endif
  
enddo


! --- CHECK CONVERGENCE ---

crystallite_todo = .false.                                                                                ! done with integration
if ( .not. (mode == 2 .and. singleRun) &                                                                  ! except for local stiffness calculation:
     .and. any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then               ! any non-local not yet converged (or broken)...
  crystallite_converged = crystallite_converged .and. crystallite_localConstitution                       ! ...restart all non-local as not converged
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 5th order Runge-Kutta Cash-Karp method with adaptive step size
! (use 5th order solution to advance = "local extrapolation")
!********************************************************************
subroutine crystallite_integrateStateRKCK45(mode,gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debugger, &
                        selectiveDebugger, &
                        verboseDebugger, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_StateLoopDistribution
use numerics, only:     rTol_crystalliteState, &
                        rTol_crystalliteTemperature, &
                        subStepSizeCryst, &
                        stepIncreaseCryst
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
                        constitutive_relevantState, &
                        constitutive_subState0, &
                        constitutive_dotState, &
                        constitutive_RKCK45dotState, &
                        constitutive_collectDotState, &
                        constitutive_dotTemperature, &
                        constitutive_microstructure

implicit none


!*** input variables ***!
integer(pInt), intent(in) ::                  mode                      ! mode of calculation; 1: central solution, 2: stiffness (by perturbation)
integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                              ii, &                     ! integration point index
                                              gg                        ! grain index

!*** output variables ***!

!*** local variables ***!
integer(pInt)                                 e, &                      ! element index in element loop
                                              i, &                      ! integration point index in ip loop
                                              g, &                      ! grain index in grain loop
                                              j, &
                                              n, &                      ! stage index in integration stage loop
                                              sizeDotState, &           ! size of dot State
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
  iIter(:,ee) = ii
  gIter(:,ee) = gg
  singleRun = .true.
else
  eIter = FEsolving_execElem(1:2)
  do e = eIter(1),eIter(2)
    iIter(:,e) = FEsolving_execIP(1:2,e)
    gIter(:,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- UPDATE DEPENDENT STATES ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                       crystallite_Fe, crystallite_Fp, g, i, e)                           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- FIRST RUNGE KUTTA STEP ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
      call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                        crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
      if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      else                                                                                                ! everything is fine
        constitutive_RKCK45dotState(1,g,i,e)%p = constitutive_dotState(g,i,e)%p                           ! initial contribution to RK slope
        RKCK45dotTemperature(1,g,i,e) = crystallite_dotTemperature(g,i,e)
      endif
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- SECOND TO SIXTH RUNGE KUTTA STEP ---

do n = 1,5


  ! --- state update ---
  
  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        sizeDotState = constitutive_sizeDotState(g,i,e)
        constitutive_dotState(g,i,e)%p = 0.0_pReal
        crystallite_dotTemperature(g,i,e) = 0.0_pReal
        do j = 1,n
          constitutive_dotState(g,i,e)%p = constitutive_dotState(g,i,e)%p + a(j,n) * constitutive_RKCK45dotState(j,g,i,e)%p
          crystallite_dotTemperature(g,i,e) = crystallite_dotTemperature(g,i,e) + a(j,n) * RKCK45dotTemperature(j,g,i,e)
        enddo
        constitutive_state(g,i,e)%p(1:sizeDotState) = constitutive_subState0(g,i,e)%p(1:sizeDotState) &
                                                      + constitutive_dotState(g,i,e)%p(1:sizeDotState) * crystallite_subdt(g,i,e)
        crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                       + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO      

  
  ! --- update dependent states ---

  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                         crystallite_Fe, crystallite_Fp, g, i, e)                         ! update dependent state variables to be consistent with basic states
      endif
   enddo; enddo; enddo
  !$OMPEND PARALLEL DO


  ! --- stress integration ---
  
  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
      if (crystallite_todo(g,i,e)) then
        if (.not. crystallite_integrateStress(mode,g,i,e,c(n))) then                                      ! fraction of original time step
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMPEND CRITICAL
          else                                                                                            ! if broken local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
          endif
        endif
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO      
  

  ! --- dot state and RK dot state---
  
  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
      if (crystallite_todo(g,i,e)) then
        call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), &
                                          crystallite_Fe, crystallite_Fp, crystallite_Temperature(g,i,e), & 
                                          c(n)*crystallite_subdt(g,i,e), g,i,e)                           ! fraction of original timestep
        crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                        crystallite_Temperature(g,i,e),g,i,e)
        if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
             .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMPEND CRITICAL
          else                                                                                            ! if broken local...
            crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
          endif
        else                                                                                              ! everything is fine
          constitutive_RKCK45dotState(n+1,g,i,e)%p = constitutive_dotState(g,i,e)%p
          RKCK45dotTemperature(n+1,g,i,e) = crystallite_dotTemperature(g,i,e)
        endif
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO

enddo  


! --- STATE UPDATE WITH ERROR ESTIMATE FOR STATE AND TEMPERATURE ---

stateResiduum = 0.0_pReal
temperatureResiduum = 0.0_pReal
relStateResiduum = 0.0_pReal
relTemperatureResiduum = 0.0_pReal
!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
      sizeDotState = constitutive_sizeDotState(g,i,e)
      
      
      ! --- absolute residuum in state and temperature ---
      
      do j = 1,6
        stateResiduum(1:sizeDotState,g,i,e) = stateResiduum(1:sizeDotState,g,i,e) &
                                        + db(j) * constitutive_RKCK45dotState(j,g,i,e)%p(1:sizeDotState) * crystallite_subdt(g,i,e)
        temperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) + db(j) * RKCK45dotTemperature(j,g,i,e) * crystallite_subdt(g,i,e)
      enddo
      
      
      ! --- dot state and dot temperature ---
      
      constitutive_dotState(g,i,e)%p = 0.0_pReal
      crystallite_dotTemperature(g,i,e) = 0.0_pReal
      do j = 1,6
        constitutive_dotState(g,i,e)%p = constitutive_dotState(g,i,e)%p + b(j) * constitutive_RKCK45dotState(j,g,i,e)%p
        crystallite_dotTemperature(g,i,e) = crystallite_dotTemperature(g,i,e) + b(j) * RKCK45dotTemperature(j,g,i,e)
      enddo


      ! --- state and temperature update and relative residui ---      

      constitutive_state(g,i,e)%p(1:sizeDotState) = constitutive_subState0(g,i,e)%p(1:sizeDotState) &
                                                    + constitutive_dotState(g,i,e)%p(1:sizeDotState) * crystallite_subdt(g,i,e)
      crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                     + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
      forall (s = 1:sizeDotState, constitutive_state(g,i,e)%p(s) > constitutive_relevantState(g,i,e)%p(s)) &
        relStateResiduum(s,g,i,e) = abs(stateResiduum(s,g,i,e)) / constitutive_state(g,i,e)%p(s) / rTol_crystalliteState
      if (crystallite_Temperature(g,i,e) > 0) &
        relTemperatureResiduum(g,i,e) = abs(temperatureResiduum(g,i,e)) / crystallite_Temperature(g,i,e) &
                                                                        / rTol_crystalliteTemperature
        
      if (verboseDebugger .and. selectiveDebugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) '::: updateState',g,i,e
          write(6,*)
          write(6,'(a,/,12(f12.1,x))') 'updateState: absolute residuum', stateResiduum(1:sizeDotState,g,i,e)
          write(6,*)
          write(6,'(a,/,12(f12.1,x))') 'updateState: resid tolerance', relStateResiduum(1:sizeDotState,g,i,e)
          write(6,*)
!          write(6,'(a)') 'updateState: RKCK45dotState'
!          do j = 1,6
!            write(6,'(12(e14.8,x))') constitutive_RKCK45dotState(j,g,i,e)%p(1:sizeDotState)
!            write(6,*)
!          enddo
          write(6,'(a,/,12(e14.8,x))') 'updateState: dotState', constitutive_dotState(g,i,e)%p(1:sizeDotState)
          write(6,*)
          write(6,'(a,/,12(e14.8,x))') 'updateState: new state', constitutive_state(g,i,e)%p(1:sizeDotState)
          write(6,*)
        !$OMPEND CRITICAL (write2out)
      endif
    
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO
      

! --- UPDATE DEPENDENT STATES IF RESIDUUM BELOW TOLERANCE ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      sizeDotState = constitutive_sizeDotState(g,i,e)
      if ( all(relStateResiduum(1:sizeDotState,g,i,e) < 1.0_pReal) .and. relTemperatureResiduum(g,i,e) < 1.0_pReal ) then        
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                         crystallite_Fe, crystallite_Fp, g, i, e)                         ! update dependent state variables to be consistent with basic states
      endif
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- FINAL STRESS INTEGRATION STEP IF RESIDUUM BELOW TOLERANCE ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
      sizeDotState = constitutive_sizeDotState(g,i,e)
      if ( all(relStateResiduum(1:sizeDotState,g,i,e) < 1.0_pReal) .and. relTemperatureResiduum(g,i,e) < 1.0_pReal ) then        

        if (crystallite_integrateStress(mode,g,i,e)) then
          crystallite_converged(g,i,e) = .true.                                                           ! ... converged per definitionem
          crystallite_todo(g,i,e) = .false.                                                               ! ... integration done
          !$OMP CRITICAL (distributionState)
            debug_StateLoopDistribution(6,mode) = debug_StateLoopDistribution(6,mode) + 1
          !$OMPEND CRITICAL (distributionState)
        else
          if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
            !$OMP CRITICAL
              crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
            !$OMPEND CRITICAL
          endif
        endif
        
      endif
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- nonlocal convergence check ---

if (verboseDebugger .and. mode==1) write(6,*) 'crystallite_converged',crystallite_converged
if ( .not. (mode == 2 .and. singleRun) ) then                                                             ! except for local stiffness calculation:
  
  if ( any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                   ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
  
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 1nd order Euler method with adaptive step size
!********************************************************************
subroutine crystallite_integrateStateAdaptiveEuler(mode,gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debugger, &
                        selectiveDebugger, &
                        verboseDebugger, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_StateLoopDistribution
use numerics, only:     rTol_crystalliteState, &
                        rTol_crystalliteTemperature, &
                        subStepSizeCryst, &
                        stepIncreaseCryst
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
                        constitutive_relevantState, &
                        constitutive_subState0, &
                        constitutive_dotState, &
                        constitutive_collectDotState, &
                        constitutive_dotTemperature, &
                        constitutive_microstructure

implicit none


!*** input variables ***!
integer(pInt), intent(in) ::                  mode                      ! mode of calculation; 1: central solution, 2: stiffness (by perturbation)
integer(pInt), optional, intent(in)::         ee, &                     ! element index
                                              ii, &                     ! integration point index
                                              gg                        ! grain index

!*** output variables ***!

!*** local variables ***!
integer(pInt)                                 e, &                      ! element index in element loop
                                              i, &                      ! integration point index in ip loop
                                              g, &                      ! grain index in grain loop
                                              j, &
                                              n, &                      ! stage index in integration stage loop
                                              sizeDotState, &           ! size of dot State
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
  iIter(:,ee) = ii
  gIter(:,ee) = gg
  singleRun = .true.
else
  eIter = FEsolving_execElem(1:2)
  do e = eIter(1),eIter(2)
    iIter(:,e) = FEsolving_execIP(1:2,e)
    gIter(:,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- UPDATE DEPENDENT STATES (EULER INTEGRATION) ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                       crystallite_Fe, crystallite_Fp, g, i, e)                           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- DOT STATE AND TEMPERATURE (EULER INTEGRATION) ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
  
      call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                        crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
  
      if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      else  
        stateResiduum(:,g,i,e) = - 0.5_pReal * constitutive_dotState(g,i,e)%p * crystallite_subdt(g,i,e)  ! contribution to absolute residuum in state and temperature
        temperatureResiduum(g,i,e) = - 0.5_pReal * crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
      endif
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- STATE UPDATE (EULER INTEGRATION) ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      sizeDotState = constitutive_sizeDotState(g,i,e)
      constitutive_state(g,i,e)%p(1:sizeDotState) = constitutive_subState0(g,i,e)%p(1:sizeDotState) &
                                                    + constitutive_dotState(g,i,e)%p(1:sizeDotState) * crystallite_subdt(g,i,e)
      crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                     + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO      


! --- UPDATE DEPENDENT STATES (EULER INTEGRATION) ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                       crystallite_Fe, crystallite_Fp, g, i, e)                         ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- STRESS INTEGRATION (EULER INTEGRATION) ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
    selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
    if (crystallite_todo(g,i,e)) then
      if (.not. crystallite_integrateStress(mode,g,i,e)) then
        if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
          !$OMPEND CRITICAL
        else                                                                                            ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO      


! --- DOT STATE AND TEMPERATURE (HEUN METHOD) ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
    selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
    if (crystallite_todo(g,i,e)) then
      call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                        crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
      if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                        ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then             ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                            ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                     ! ...all non-locals skipped
          !$OMPEND CRITICAL
        else                                                                                            ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                             ! ... skip this one next time
        endif
      endif      
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- ERROR ESTIMATE FOR STATE AND TEMPERATURE (HEUN METHOD) ---

relStateResiduum = 0.0_pReal
relTemperatureResiduum = 0.0_pReal
!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode==1)
      sizeDotState = constitutive_sizeDotState(g,i,e)      
      
      
      ! --- contribution of heun step to absolute residui ---
      
      stateResiduum(:,g,i,e) = stateResiduum(:,g,i,e) &
                               + 0.5_pReal * constitutive_dotState(g,i,e)%p * crystallite_subdt(g,i,e)    ! contribution to absolute residuum in state and temperature
      temperatureResiduum(g,i,e) = temperatureResiduum(g,i,e) &
                                   + 0.5_pReal * crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)


      ! --- relative residui ---      

      forall (s = 1:sizeDotState, constitutive_state(g,i,e)%p(s) > constitutive_relevantState(g,i,e)%p(s)) &
        relStateResiduum(s,g,i,e) = abs(stateResiduum(s,g,i,e)) / constitutive_state(g,i,e)%p(s) / rTol_crystalliteState
      if (crystallite_Temperature(g,i,e) > 0) &
        relTemperatureResiduum(g,i,e) = abs(temperatureResiduum(g,i,e)) / crystallite_Temperature(g,i,e) &
                                                                        / rTol_crystalliteTemperature
        
      if (verboseDebugger .and. selectiveDebugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) '::: updateState',g,i,e
          write(6,*)
          write(6,'(a,/,12(f12.1,x))') 'updateState: absolute residuum', stateResiduum(1:sizeDotState,g,i,e)
          write(6,*)
          write(6,'(a,/,12(f12.1,x))') 'updateState: resid tolerance', relStateResiduum(1:sizeDotState,g,i,e)
          write(6,*)
          write(6,'(a,/,12(e14.8,x))') 'updateState: dotState', constitutive_dotState(g,i,e)%p(1:sizeDotState) &
                                                      - 2.0_pReal * stateResiduum(1:sizeDotState,g,i,e) / crystallite_subdt(g,i,e)  ! calculate former dotstate from higher order solution and state residuum
          write(6,*)
          write(6,'(a,/,12(e14.8,x))') 'updateState: new state', constitutive_state(g,i,e)%p(1:sizeDotState)
          write(6,*)
        !$OMPEND CRITICAL (write2out)
      endif
            
      
      ! --- converged ? ---

      if ( all(relStateResiduum(1:sizeDotState,g,i,e) < 1.0_pReal) .and. relTemperatureResiduum(g,i,e) < 1.0_pReal ) then        

        crystallite_converged(g,i,e) = .true.                                                           ! ... converged per definitionem
        crystallite_todo(g,i,e) = .false.                                                               ! ... integration done
        !$OMP CRITICAL (distributionState)
          debug_StateLoopDistribution(6,mode) = debug_StateLoopDistribution(6,mode) + 1
        !$OMPEND CRITICAL (distributionState)

      endif

    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- NONLOCAL CONVERGENCE CHECK ---

if (verboseDebugger .and. mode==1) write(6,*) 'crystallite_converged',crystallite_converged
if ( .not. (mode == 2 .and. singleRun) ) then                                                             ! except for local stiffness calculation:  
  if ( any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then                   ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! 1st order explicit Euler method 
!********************************************************************
subroutine crystallite_integrateStateEuler(mode,gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debugger, &
                        selectiveDebugger, &
                        verboseDebugger, &
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
integer(pInt), intent(in) ::                  mode                      ! mode of calculation; 1: central solution, 2: stiffness (by perturbation)
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
logical                                       singleRun                 ! flag indicating computation for single (g,i,e) triple


if (present(ee) .and. present(ii) .and. present(gg)) then
  eIter = ee
  iIter(:,ee) = ii
  gIter(:,ee) = gg
  singleRun = .true.
else
  eIter = FEsolving_execElem(1:2)
  do e = eIter(1),eIter(2)
    iIter(:,e) = FEsolving_execIP(1:2,e)
    gIter(:,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --- UPDATE DEPENDENT STATES ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                       crystallite_Fe, crystallite_Fp, g, i, e)                           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- DOT STATE AND TEMPERATURE ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
      call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                        crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g,i,e)
      crystallite_dotTemperature(g,i,e) = constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e), &
                                                                      crystallite_Temperature(g,i,e),g,i,e)
      if ( any(constitutive_dotState(g,i,e)%p/=constitutive_dotState(g,i,e)%p) &                          ! NaN occured in dotState
           .or. crystallite_dotTemperature(g,i,e)/=crystallite_dotTemperature(g,i,e) ) then               ! NaN occured in dotTemperature
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        else                                                                                              ! if broken local...
          crystallite_todo(g,i,e) = .false.                                                               ! ... skip this one next time
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- UPDATE STATE AND TEMPERATURE ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
      mySizeDotState = constitutive_sizeDotState(g,i,e)
      constitutive_state(g,i,e)%p(1:mySizeDotState) = constitutive_subState0(g,i,e)%p(1:mySizeDotState) &
                                                    + constitutive_dotState(g,i,e)%p(1:mySizeDotState) * crystallite_subdt(g,i,e)
      crystallite_Temperature(g,i,e) = crystallite_subTemperature0(g,i,e) &
                                        + crystallite_dotTemperature(g,i,e) * crystallite_subdt(g,i,e)
                                        
      if (verboseDebugger .and. selectiveDebugger) then
        !$OMP CRITICAL (write2out)
          write(6,*) '::: updateState',g,i,e
          write(6,*)
          write(6,'(a,/,12(e14.8,x))') 'updateState: dotState', constitutive_dotState(g,i,e)%p(1:mySizeDotState)
          write(6,*)
          write(6,'(a,/,12(e14.8,x))') 'updateState: new state', constitutive_state(g,i,e)%p(1:mySizeDotState)
          write(6,*)
        !$OMPEND CRITICAL (write2out)
      endif
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO      

  
! --- UPDATE DEPENDENT STATES ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                       crystallite_Fe, crystallite_Fp, g, i, e)                           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- STRESS INTEGRATION ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
    if (crystallite_todo(g,i,e)) then
      if (crystallite_integrateStress(mode,g,i,e)) then
        crystallite_converged(g,i,e) = .true.
        !$OMP CRITICAL (distributionState)
          debug_StateLoopDistribution(1,mode) = debug_StateLoopDistribution(1,mode) + 1
        !$OMPEND CRITICAL (distributionState)
      else                                                                                                ! broken stress integration
        if (.not. crystallite_localConstitution(g,i,e)) then                                              ! if broken non-local...
          !$OMP CRITICAL
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        endif
      endif
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO      


! --- CHECK NON-LOCAL CONVERGENCE ---

crystallite_todo = .false.                                                                                ! done with integration
if ( .not. (mode == 2 .and. singleRun) &                                                                  ! except for local stiffness calculation:
     .and. any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then               ! any non-local not yet converged (or broken)...
  crystallite_converged = crystallite_converged .and. crystallite_localConstitution                       ! ...restart all non-local as not converged
endif

endsubroutine



!********************************************************************
! integrate stress, state and Temperature with 
! adaptive 1st order explicit Euler method 
! using Fixed Point Iteration to adapt the stepsize  
!********************************************************************
subroutine crystallite_integrateStateFPI(mode,gg,ii,ee)

!*** variables and functions from other modules ***!
use prec, only:         pInt, &
                        pReal
use debug, only:        debugger, &
                        selectiveDebugger, &
                        verboseDebugger, &
                        debug_e, &
                        debug_i, &
                        debug_g, &
                        debug_StateLoopDistribution
use numerics, only:     nState
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
integer(pInt), intent(in) ::                  mode                      ! mode of calculation; 1: central solution, 2: stiffness (by perturbation)
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
logical                                       singleRun                 ! flag indicating computation for single (g,i,e) triple


if (present(ee) .and. present(ii) .and. present(gg)) then
  eIter = ee
  iIter(:,ee) = ii
  gIter(:,ee) = gg
  singleRun = .true.
else
  eIter = FEsolving_execElem(1:2)
  do e = eIter(1),eIter(2)
    iIter(:,e) = FEsolving_execIP(1:2,e)
    gIter(:,e) = (/1,homogenization_Ngrains(mesh_element(3,e))/)
  enddo
  singleRun = .false.
endif


! --+>> PREGUESS FOR STATE <<+--

! --- UPDATE DEPENDENT STATES ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    if (crystallite_todo(g,i,e)) then
      call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                       crystallite_Fe, crystallite_Fp, g, i, e)                           ! update dependent state variables to be consistent with basic states
    endif
 enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- DOT STATES ---

!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
    if (crystallite_todo(g,i,e)) then
      constitutive_previousDotState2(g,i,e)%p = 0.0_pReal
      constitutive_previousDotState(g,i,e)%p = 0.0_pReal
      call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                        crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g, i, e)
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --- STATE & TEMPERATURE UPDATE ---

crystallite_statedamper = 1.0_pReal
!$OMP PARALLEL DO
  do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                          ! iterate over elements, ips and grains
    selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
    if (crystallite_todo(g,i,e)) then
      crystallite_stateConverged(g,i,e) = crystallite_updateState(g,i,e)                                  ! update state
      crystallite_temperatureConverged(g,i,e) = crystallite_updateTemperature(g,i,e)                      ! update temperature
      if ( .not. crystallite_localConstitution(g,i,e) .and. .not. crystallite_todo(g,i,e) ) then          ! if updateState or updateTemperature signals broken non-local... 
        !$OMP CRITICAL 
          crystallite_todo = crystallite_todo .and. crystallite_localConstitution                         ! ...all non-locals skipped
        !$OMPEND CRITICAL
      endif
    endif
  enddo; enddo; enddo
!$OMPEND PARALLEL DO


! --+>> STATE LOOP <<+--

NiterationState = 0_pInt

do while (any(crystallite_todo) .and. NiterationState < nState )                                          ! convergence loop for crystallite
 
  NiterationState = NiterationState + 1_pInt


  ! --- STRESS INTEGRATION ---

  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
      if (crystallite_todo(g,i,e)) then
        crystallite_todo(g,i,e) = crystallite_integrateStress(mode,g,i,e)
        if ( .not. crystallite_localConstitution(g,i,e) .and. .not. crystallite_todo(g,i,e)) then         ! if broken non-local... 
          !$OMP CRITICAL 
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        endif
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO

  if (verboseDebugger .and. mode == 1) then
    !$OMP CRITICAL (write2out)
      write(6,*) count(crystallite_todo(:,:,:)),'grains todo after stress integration'
    !$OMPEND CRITICAL (write2out)
  endif


  ! --- DOT STATES ---
  
  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
      if (crystallite_todo(g,i,e)) then
        constitutive_previousDotState2(g,i,e)%p = constitutive_previousDotState(g,i,e)%p                  ! wind forward dotStates
        constitutive_previousDotState(g,i,e)%p = constitutive_dotState(g,i,e)%p
        call constitutive_collectDotState(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, &
                                          crystallite_Fp, crystallite_Temperature(g,i,e), crystallite_subdt(g,i,e), g, i, e)
      endif
  enddo; enddo; enddo
  !$OMPEND PARALLEL DO


  ! --- STATE & TEMPERATURE UPDATE ---

  crystallite_statedamper = 1.0_pReal
  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      selectiveDebugger = (e == debug_e .and. i == debug_i .and. g == debug_g .and. mode == 1)
      if (crystallite_todo(g,i,e)) then
        
        ! --- state damper ---
        
        dot_prod12 = dot_product( constitutive_dotState(g,i,e)%p - constitutive_previousDotState(g,i,e)%p, &
                              constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p )
        dot_prod22 = dot_product( constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p, &
                                  constitutive_previousDotState(g,i,e)%p - constitutive_previousDotState2(g,i,e)%p )
        if (      dot_prod22 > 0.0_pReal &
            .and. (     dot_prod12 < 0.0_pReal &
                   .or. dot_product(constitutive_dotState(g,i,e)%p, constitutive_previousDotState(g,i,e)%p) < 0.0_pReal) ) &
          crystallite_statedamper(g,i,e) = 0.75_pReal + 0.25_pReal * tanh(2.0_pReal + 4.0_pReal * dot_prod12 / dot_prod22)

        ! --- updates ---
        
        crystallite_stateConverged(g,i,e) = crystallite_updateState(g,i,e)                                ! update state
        crystallite_temperatureConverged(g,i,e) = crystallite_updateTemperature(g,i,e)                    ! update temperature
        crystallite_converged(g,i,e) = crystallite_stateConverged(g,i,e) .and. crystallite_temperatureConverged(g,i,e)
        if ( .not. crystallite_localConstitution(g,i,e) .and. .not. crystallite_todo(g,i,e)) then         ! if updateState or updateTemperature signals broken non-local... 
          !$OMP CRITICAL 
            crystallite_todo = crystallite_todo .and. crystallite_localConstitution                       ! ...all non-locals skipped
          !$OMPEND CRITICAL
        elseif (crystallite_converged(g,i,e)) then
          !$OMP CRITICAL (distributionState)
            debug_StateLoopDistribution(NiterationState,mode) = debug_StateLoopDistribution(NiterationState,mode) + 1
          !$OMPEND CRITICAL (distributionState)
        endif
      endif
    enddo; enddo; enddo
  !$OMPEND PARALLEL DO


  ! --- UPDATE DEPENDENT STATES ---

  !$OMP PARALLEL DO
    do e=eIter(1),eIter(2); do i=iIter(1,e),iIter(2,e); do g=gIter(1,e),gIter(2,e)                        ! iterate over elements, ips and grains
      if (crystallite_todo(g,i,e)) then
        call constitutive_microstructure(crystallite_Temperature(g,i,e), crystallite_Tstar_v(:,g,i,e), &
                                         crystallite_Fe, crystallite_Fp, g, i, e)                         ! update dependent state variables to be consistent with basic states
      endif
   enddo; enddo; enddo
  !$OMPEND PARALLEL DO
  
  if (verboseDebugger .and. mode == 1) then
    !$OMP CRITICAL (write2out)
      write(6,*) count(crystallite_converged(:,:,:)),'grains converged after state integration no.', NiterationState
      write(6,*)
    !$OMPEND CRITICAL (write2out)
  endif

  
  ! --- CONVERGENCE CHECK ---

  if ( .not. (mode == 2 .and. singleRun) &                                                                ! except for local stiffness calculation:
       .and. any(.not. crystallite_converged .and. .not. crystallite_localConstitution)) then             ! any non-local not yet converged (or broken)...
    crystallite_converged = crystallite_converged .and. crystallite_localConstitution                     ! ...restart all non-local as not converged
  endif
  
  crystallite_todo = crystallite_todo .and. .not. crystallite_converged                                   ! skip all converged
  
  if (verboseDebugger .and. mode == 1) then
    !$OMP CRITICAL (write2out)
      write(6,*) count(crystallite_converged(:,:,:)),'grains converged after non-local check'
      write(6,*) count(crystallite_todo(:,:,:)),'grains todo after state integration no.', NiterationState
      write(6,*)
    !$OMPEND CRITICAL (write2out)
  endif
  
enddo                                                                                           ! crystallite convergence loop  

endsubroutine



!********************************************************************
! update the internal state of the constitutive law
! and tell whether state has converged
!********************************************************************
function crystallite_updateState(g,i,e)

!*** variables and functions from other modules ***!
use prec, only:           pReal, &
                          pInt, &
                          pLongInt
use numerics, only:       rTol_crystalliteState
use constitutive, only:   constitutive_dotState, &
                          constitutive_previousDotState, &
                          constitutive_sizeDotState, &
                          constitutive_subState0, &
                          constitutive_state, &
                          constitutive_relevantState, &
                          constitutive_microstructure
use debug, only:          debugger, &
                          selectiveDebugger, &
                          verboseDebugger

!*** input variables ***!
integer(pInt), intent(in):: e, &                     ! element index
                            i, &                     ! integration point index
                            g                        ! grain index

!*** output variables ***!
logical                     crystallite_updateState  ! flag indicating if integration suceeded

!*** local variables ***!
real(pReal), dimension(constitutive_sizeDotState(g,i,e)) :: residuum ! residuum from evolution of microstructure
integer(pInt)               mySize


mySize = constitutive_sizeDotState(g,i,e)

! correct my dotState 
constitutive_dotState(g,i,e)%p(1:mySize) = constitutive_dotState(g,i,e)%p(1:mySize) * crystallite_statedamper(g,i,e) &
                                 + constitutive_previousDotState(g,i,e)%p(1:mySize) * (1.0_pReal-crystallite_statedamper(g,i,e))

residuum = constitutive_state(g,i,e)%p(1:mySize) - constitutive_subState0(g,i,e)%p(1:mySize) &
                                                - constitutive_dotState(g,i,e)%p(1:mySize) * crystallite_subdt(g,i,e)

if (any(residuum/=residuum)) then                                    ! if NaN occured then return without changing the state...
  crystallite_updateState = .false.                                  ! ...indicate state update failed
  crystallite_todo(g,i,e) = .false.                                  ! ...no need to calculate any further
  if (verboseDebugger) then
    !$OMP CRITICAL (write2out)
      write(6,*) '::: updateState encountered NaN',g,i,e
    !$OMPEND CRITICAL (write2out)
  endif   
  return
endif

constitutive_state(g,i,e)%p(1:mySize) = constitutive_state(g,i,e)%p(1:mySize) - residuum

! setting flag to true if state is below relative tolerance, otherwise set it to false
crystallite_updateState = all(    constitutive_state(g,i,e)%p(1:mySize) < constitutive_relevantState(g,i,e)%p(1:mySize) &
                              .or. abs(residuum) < rTol_crystalliteState*abs(constitutive_state(g,i,e)%p(1:mySize)) )

if (verboseDebugger .and. selectiveDebugger) then
  !$OMP CRITICAL (write2out)
    if (crystallite_updateState) then
      write(6,*) '::: updateState converged',g,i,e
    else
      write(6,*) '::: updateState did not converge',g,i,e
    endif
    write(6,*)
    write(6,'(a,f6.1)') 'updateState: crystallite_statedamper',crystallite_statedamper(g,i,e)
    write(6,*)
    write(6,'(a,/,12(e14.8,x))') 'updateState: dotState',constitutive_dotState(g,i,e)%p(1:mySize)
    write(6,*)
    write(6,'(a,/,12(e14.8,x))') 'updateState: new state',constitutive_state(g,i,e)%p(1:mySize)
    write(6,*)
    write(6,'(a,/,12(f12.1,x))') 'updateState: resid tolerance',abs(residuum / rTol_crystalliteState &
                                                                             / constitutive_state(g,i,e)%p(1:mySize))
    write(6,*)
  !$OMPEND CRITICAL (write2out)
endif

endfunction



!********************************************************************
! update the temperature of the grain
! and tell whether it has converged
!********************************************************************
 function crystallite_updateTemperature(&
   g,&              ! grain number
   i,&              ! integration point number
   e &              ! element number
 )
 
 !*** variables and functions from other modules ***!
 use prec, only:                      pReal, &
                                      pInt, &
                                      pLongInt
 use numerics, only:                  rTol_crystalliteTemperature
 use constitutive, only:              constitutive_dotTemperature
 use debug, only:                     debugger
 
 !*** input variables ***!
 integer(pInt), intent(in)::          e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 
 !*** output variables ***!
 logical                              crystallite_updateTemperature ! flag indicating if integration suceeded

 !*** local variables ***!
 real(pReal) residuum                                               ! residuum from evolution of temperature
 
 ! calculate the residuum 
 residuum = crystallite_Temperature(g,i,e) - crystallite_subTemperature0(g,i,e) - &
            crystallite_subdt(g,i,e) * &
            constitutive_dotTemperature(crystallite_Tstar_v(:,g,i,e),crystallite_Temperature(g,i,e),g,i,e)
 
 ! if NaN occured then return without changing the state
 if (residuum/=residuum) then
   crystallite_updateTemperature = .false.                                  ! indicate update failed
   crystallite_todo(g,i,e) = .false.                                        ! ...no need to calculate any further
   !$OMP CRITICAL (write2out)
   write(6,*) '::: updateTemperature encountered NaN',g,i,e
   !$OMPEND CRITICAL (write2out)
   return
 endif
 
 ! update the microstructure
 crystallite_Temperature(g,i,e) = crystallite_Temperature(g,i,e) - residuum
 
 ! setting flag to true if residuum is below relative tolerance (or zero Kelvin), otherwise set it to false
 crystallite_updateTemperature = crystallite_Temperature(g,i,e) == 0.0_pReal .or. &
                                  abs(residuum) < rTol_crystalliteTemperature*crystallite_Temperature(g,i,e)
 
 return

 endfunction



!***********************************************************************
!***     calculation of stress (P) with time integration             ***
!***     based on a residuum in Lp and intermediate                  ***
!***     acceleration of the Newton-Raphson correction               ***
!***********************************************************************
 function crystallite_integrateStress(&
     mode, &      ! 1: central solution, 2: stiffness (by perturbation)
     g,&          ! grain number
     i,&          ! integration point number
     e,&          ! element number
     fraction &
     )
     

 !*** variables and functions from other modules ***!
 use prec, only:                      pReal, &
                                      pInt, &
                                      pLongInt
 use numerics, only:                  nStress, &
                                      aTol_crystalliteStress, &
                                      rTol_crystalliteStress, &
                                      iJacoLpresiduum, &
                                      relevantStrain
 use debug, only:                     debugger, &
                                      selectiveDebugger, &
                                      verboseDebugger, &
                                      debug_cumLpCalls, &
                                      debug_cumLpTicks, &
                                      debug_StressLoopDistribution, &
                                      debug_LeapfrogBreakDistribution
 use constitutive, only:              constitutive_homogenizedC, &
                                      constitutive_LpAndItsTangent
 use math, only:                      math_mul33x33, &
                                      math_mul66x6, &
                                      math_mul99x99, &
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
 integer(pInt), intent(in)::          mode, &                       ! 1 or 2
                                      e, &                          ! element index
                                      i, &                          ! integration point index
                                      g                             ! grain index
 real(pReal), optional, intent(in) :: fraction                      ! fraction of timestep

 !*** output variables ***!
 logical                              crystallite_integrateStress   ! flag indicating if integration suceeded
 
 !*** local variables ***!
 real(pReal), dimension(3,3)::        Fg_new, &                     ! deformation gradient at end of timestep
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
 real(pReal), dimension(6)::          Tstar_v                       ! 2nd Piola-Kirchhoff Stress in Mandel-Notation
 real(pReal), dimension(9,9)::        dLpdT_constitutive, &         ! partial derivative of plastic velocity gradient calculated by constitutive law
                                      dTdLp, &                      ! partial derivative of 2nd Piola-Kirchhoff stress
                                      dRdLp, &                      ! partial derivative of residuum (Jacobian for NEwton-Raphson scheme)
                                      invdRdLp                      ! inverse of dRdLp
 real(pReal), dimension(3,3,3,3)::    C                             ! 4th rank elasticity tensor
 real(pReal), dimension(6,6)::        C_66                          ! simplified 2nd rank elasticity tensor 
 real(pReal)                          p_hydro, &                    ! volumetric part of 2nd Piola-Kirchhoff Stress
                                      det, &                        ! determinant
                                      leapfrog, &                   ! acceleration factor for Newton-Raphson scheme
                                      maxleap, &                    ! maximum acceleration factor
                                      dt                            ! time increment
 logical                              error                         ! flag indicating an error
 integer(pInt)                        NiterationStress, &           ! number of stress integrations
                                      dummy, &
                                      h, &
                                      j, &
                                      k, &
                                      l, &
                                      m, &
                                      n, &
                                      jacoCounter                   ! counter to check for Jacobian update
 integer(pLongInt)                    tick, &
                                      tock, &
                                      tickrate, &
                                      maxticks
 
 ! be pessimistic
 crystallite_integrateStress = .false.

 ! only integrate over fraction of timestep?
 if (present(fraction)) then
   dt = crystallite_subdt(g,i,e) * fraction
   Fg_new = crystallite_subF0(:,:,g,i,e) + (crystallite_subF(:,:,g,i,e) - crystallite_subF0(:,:,g,i,e)) * fraction
 else     
   dt = crystallite_subdt(g,i,e)
   Fg_new = crystallite_subF(:,:,g,i,e)
 endif

 ! feed local variables
 Fp_current =   crystallite_subFp0(:,:,g,i,e)
 Tstar_v =      crystallite_Tstar_v(:,g,i,e)
 Lpguess_old =  crystallite_Lp(:,:,g,i,e)                           ! consider present Lp good (i.e. worth remembering) ...
 Lpguess =      crystallite_Lp(:,:,g,i,e)                           ! ... and take it as first guess

 
 ! inversion of Fp_current...
 invFp_current = math_inv3x3(Fp_current)                            
 if (all(invFp_current == 0.0_pReal)) then                          ! ... failed?
   if (verboseDebugger .and. selectiveDebugger) then 
     !$OMP CRITICAL (write2out)
       write(6,*) '::: integrateStress failed on invFp_current inversion',g,i,e
       write(6,*)
       write(6,'(a11,i3,x,i2,x,i5,/,3(3(f12.7,x)/))') 'invFp_new at ',g,i,e,invFp_new
     !$OMPEND CRITICAL (write2out)
   endif
   return
 endif
 
 A = math_mul33x33(transpose(invFp_current), math_mul33x33(transpose(Fg_new),math_mul33x33(Fg_new,invFp_current)))
 
 ! get elasticity tensor
 C_66 = constitutive_homogenizedC(g,i,e)
! if (debugger) write(6,'(a,/,6(6(f10.4,x)/))') 'elasticity',C_66(1:6,:)/1e9
 C = math_Mandel66to3333(C_66)
 
 ! start LpLoop with no acceleration
 NiterationStress = 0_pInt
 leapfrog = 1.0_pReal
 maxleap = 1024.0_pReal
 jacoCounter = 0_pInt

LpLoop: do
   
   ! increase loop counter
   NiterationStress = NiterationStress + 1
   
   ! too many loops required ?
   if (NiterationStress > nStress) then
     if (verboseDebugger) then 
       !$OMP CRITICAL (write2out)
         write(6,*) '::: integrateStress reached loop limit at ',g,i,e
         write(6,*)
       !$OMPEND CRITICAL (write2out)
     endif
     return
   endif
   
   B = math_I3 - dt*Lpguess
   BT = transpose(B)
   AB = math_mul33x33(A,B)
   BTA = math_mul33x33(BT,A)
   
   ! calculate 2nd Piola-Kirchhoff stress tensor
   Tstar_v = 0.5_pReal*math_mul66x6(C_66,math_mandel33to6(math_mul33x33(BT,AB)-math_I3))
   p_hydro = sum(Tstar_v(1:3))/3.0_pReal
   forall(n=1:3) Tstar_v(n) = Tstar_v(n) - p_hydro                  ! get deviatoric stress tensor
   
   ! calculate plastic velocity gradient and its tangent according to constitutive law
   call system_clock(count=tick,count_rate=tickrate,count_max=maxticks)
   call constitutive_LpAndItsTangent(Lp_constitutive, dLpdT_constitutive, Tstar_v, crystallite_Temperature(g,i,e), g, i, e)
   call system_clock(count=tock,count_rate=tickrate,count_max=maxticks)
   debug_cumLpCalls = debug_cumLpCalls + 1_pInt
   debug_cumLpTicks  = debug_cumLpTicks + tock-tick
   if (tock < tick) debug_cumLpTicks = debug_cumLpTicks + maxticks
   if (verboseDebugger .and. selectiveDebugger) then
     !$OMP CRITICAL (write2out)
       write(6,'(a,i3,x,i2,x,i5,x,a,x,i3)') '::: integrateStress at ' ,g,i,e, ' ; iteration ', NiterationStress
       write(6,*)
       write(6,'(a,/,3(3(e20.7,x)/))') 'Lp_constitutive', Lp_constitutive
       write(6,'(a,/,3(3(e20.7,x)/))') 'Lpguess', Lpguess
     !$OMPEND CRITICAL (write2out)
   endif

   ! update current residuum
   residuum = Lpguess - Lp_constitutive

   ! Check for convergence of loop
   if (.not.(any(residuum/=residuum)) .and. &                       ! exclude any NaN in residuum
       ( maxval(abs(residuum)) < aTol_crystalliteStress .or. &      ! below absolute tolerance .or.
         ( any(abs(dt*Lpguess) > relevantStrain) .and. &            ! worth checking? .and.
             maxval(abs(residuum/Lpguess), abs(dt*Lpguess) > relevantStrain) < rTol_crystalliteStress & ! below relative tolerance
         ) &
       ) &
      ) &
     exit LpLoop
   
   ! NaN occured at regular speed?
   if (any(residuum/=residuum) .and. leapfrog == 1.0) then
     if (debugger) then 
       !$OMP CRITICAL (write2out)
         write(6,'(a,i3,x,i2,x,i5,x,a,i3,x,a)') '::: integrateStress encountered NaN at ',g,i,e,&
                                                '; iteration ', NiterationStress, &
                                                '>> returning..!'
       !$OMPEND CRITICAL (write2out)
     endif
     return

   ! something went wrong at accelerated speed?
   elseif (leapfrog > 1.0_pReal .and. &                             ! at fast pace .and.
            ( sum(residuum*residuum) > sum(residuum_old*residuum_old) .or. &  ! worse residuum .or.
              sum(residuum*residuum_old) < 0.0_pReal .or. &         ! residuum changed sign (overshoot) .or.
              any(residuum/=residuum) &                             ! NaN occured
            ) &
          ) then
      if (verboseDebugger) then 
       !$OMP CRITICAL (write2out)
         write(6,'(a,i3,x,i2,x,i5,x,a,i3)') '::: integrateStress encountered high-speed crash at ',g,i,e,&
                                            '; iteration ', NiterationStress
       !$OMPEND CRITICAL (write2out)
     endif
     maxleap = 0.5_pReal * leapfrog                                 ! limit next acceleration
     leapfrog = 1.0_pReal                                           ! grinding halt
     jacoCounter = 0_pInt                                           ! reset counter for Jacobian update (we want to do an update next time!)
     
     ! restore old residuum and Lp
     Lpguess = Lpguess_old                                       
     residuum  = residuum_old
     
     debug_LeapfrogBreakDistribution(NiterationStress,mode) = debug_LeapfrogBreakDistribution(NiterationStress,mode) + 1
     
   ! residuum got better
   else
     ! calculate Jacobian for correction term 
     if (mod(jacoCounter, iJacoLpresiduum) == 0_pInt) then
       dTdLp = 0.0_pReal
       do h=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3
!       forall (h=1:3,j=1:3,k=1:3,l=1:3,m=1:3) &
         dTdLp(3*(h-1)+j,3*(k-1)+l) = dTdLp(3*(h-1)+j,3*(k-1)+l) + C(h,j,l,m)*AB(k,m)+C(h,j,m,l)*BTA(m,k)
       enddo; enddo; enddo; enddo; enddo
       dTdLp = -0.5_pReal*dt*dTdLp
       dRdLp = math_identity2nd(9) - math_mul99x99(dLpdT_constitutive,dTdLp)
       invdRdLp = 0.0_pReal
       call math_invert(9,dRdLp,invdRdLp,dummy,error)               ! invert dR/dLp --> dLp/dR
       if (error) then
         if (verboseDebugger .and. selectiveDebugger) then
           !$OMP CRITICAL (write2out)
             write(6,'(a,i3,x,i2,x,i5,x,a,i3)') '::: integrateStress failed on dR/dLp inversion at ',g,i,e, &
                                                '; iteration ', NiterationStress
             write(6,*)
             write(6,'(a,/,9(9(e15.3,x)/))') 'dRdLp',dRdLp
             write(6,'(a,/,9(9(e15.3,x)/))') 'dLpdT_constitutive',dLpdT_constitutive
             write(6,'(a,/,3(3(e20.7,x)/))') 'Lp_constitutive',Lp_constitutive
             write(6,'(a,/,3(3(e20.7,x)/))') 'Lpguess',Lpguess
           !$OMPEND CRITICAL (write2out)
         endif
         return
       else
         if (verboseDebugger .and. selectiveDebugger) then
           !$OMP CRITICAL (write2out)
             write(6,'(a,i3,x,i2,x,i5,x,a,i3)') '::: integrateStress did dR/dLp inversion at ',g,i,e, &
                                                '; iteration ', NiterationStress
             write(6,*)
             write(6,'(a,/,9(9(e15.3,x)/))') 'dRdLp',dRdLp
             write(6,'(a,/,9(9(e15.3,x)/))') 'dLpdT_constitutive',dLpdT_constitutive
           !$OMPEND CRITICAL (write2out)
         endif
       endif
     endif
     jacoCounter = jacoCounter + 1_pInt                             ! increase counter for jaco update
     
     ! remember current residuum and Lpguess
     residuum_old = residuum
     Lpguess_old = Lpguess 
     
     ! accelerate?
     if (NiterationStress > 1 .and. leapfrog < maxleap) leapfrog = 2.0_pReal * leapfrog
   endif

   ! leapfrog to updated Lp
   do k=1,3; do l=1,3; do m=1,3; do n=1,3
     Lpguess(k,l) = Lpguess(k,l) - leapfrog*invdRdLp(3*(k-1)+l,3*(m-1)+n)*residuum(m,n)
   enddo; enddo; enddo; enddo
 enddo LpLoop

 ! calculate new plastic and elastic deformation gradient
 invFp_new = math_mul33x33(invFp_current,B)
 invFp_new = invFp_new/math_det3x3(invFp_new)**(1.0_pReal/3.0_pReal)  ! regularize by det
 call math_invert3x3(invFp_new,Fp_new,det,error)
 if (error) then
   if (verboseDebugger .and. selectiveDebugger) then
     !$OMP CRITICAL (write2out)
       write(6,'(a,i3,x,i2,x,i5,x,a,x,i3)') '::: integrateStress failed on invFp_new inversion at ',g,i,e, &
                                            ' ; iteration ', NiterationStress
       write(6,*)
       write(6,'(a11,3(i3,x),/,3(3(f12.7,x)/))') 'invFp_new at ',g,i,e,invFp_new
     !$OMPEND CRITICAL (write2out)
   endif
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
 crystallite_invFp(:,:,g,i,e) = invFp_new

 ! set return flag to true
 crystallite_integrateStress = .true.
 if (verboseDebugger .and. selectiveDebugger) then 
   !$OMP CRITICAL (write2out)
   write(6,'(a,i3,x,i2,x,i5,x,a,x,i3)') '::: integrateStress converged at ',g,i,e,' ; iteration ', NiterationStress
   write(6,*)
   write(6,'(a,/,3(3(f12.7,x)/))') 'P / MPa',crystallite_P(:,:,g,i,e)/1e6
   write(6,'(a,/,3(3(f12.7,x)/))') 'Cauchy / MPa',math_mul33x33(crystallite_P(:,:,g,i,e),transpose(Fg_new))/1e6/math_det3x3(Fg_new)
   write(6,'(a,/,3(3(f12.7,x)/))') 'Fe Lp Fe^-1',math_mul33x33(Fe_new,math_mul33x33(crystallite_Lp(:,:,g,i,e),math_inv3x3(Fe_new)))
   write(6,'(a,/,3(3(f12.7,x)/))') 'Fp',crystallite_Fp(:,:,g,i,e)
   !$OMP END CRITICAL (write2out)
 endif

 !$OMP CRITICAL (distributionStress)
 debug_StressLoopDistribution(NiterationStress,mode) = debug_StressLoopDistribution(NiterationStress,mode) + 1
 !$OMPEND CRITICAL (distributionStress)

 return

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
use debug, only:                      debugger, &
                                      debug_e, debug_i, debug_g, &
                                      verboseDebugger, &
                                      selectiveDebugger
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
logical error

! --- CALCULATE ORIENTATION AND LATTICE ROTATION ---

!$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      do g = 1,homogenization_Ngrains(mesh_element(3,e))
        
        call math_pDecomposition(crystallite_Fe(:,:,g,i,e), U, R, error)                                  ! polar decomposition of Fe
        if (error) then
          call IO_warning(650, e, i, g)
          crystallite_orientation(:,g,i,e) = (/1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal/)               ! fake orientation
        else
          crystallite_orientation(:,g,i,e) = math_RtoQuaternion(transpose(R))
        endif
        
        crystallite_rotation(:,g,i,e) = &
          math_QuaternionDisorientation( math_qConj(crystallite_orientation(:,g,i,e)), &                  ! calculate grainrotation
                                         math_qConj(crystallite_orientation0(:,g,i,e)), &
                                         0_pInt )                                                         ! we don't want symmetry here  
        
      enddo
    enddo
  enddo
!$OMPEND PARALLEL DO


! --- UPDATE SOME ADDITIONAL VARIABLES THAT ARE NEEDED FOR NONLOCAL MATERIAL ---
! --- we use crystallite_orientation from above, so need a seperate loop

!$OMP PARALLEL DO
  do e = FEsolving_execElem(1),FEsolving_execElem(2)
    do i = FEsolving_execIP(1,e),FEsolving_execIP(2,e)
      selectiveDebugger = (e == debug_e .and. i == debug_i)
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
                  math_QuaternionDisorientation( crystallite_orientation(:,1,i,e), &
                                                 crystallite_orientation(:,1,neighboring_i,neighboring_e), & 
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
!$OMPEND PARALLEL DO

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
                                      math_I3, &
                                      inDeg, &
                                      math_Mandel6to33
 use mesh, only:                      mesh_element
 use material, only:                  microstructure_crystallite, &
                                      crystallite_Noutput, &
                                      material_phase, &
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
 integer(pInt)                        k,l,o,c,crystID,mySize
 logical                              error

 crystID = microstructure_crystallite(mesh_element(4,e))

 crystallite_postResults = 0.0_pReal
 c = 0_pInt
 crystallite_postResults(c+1) = crystallite_sizePostResults(crystID); c = c+1_pInt         ! size of results from cryst
 
 do o = 1,crystallite_Noutput(crystID)
   select case(crystallite_output(o,crystID))
     case ('phase')
       crystallite_postResults(c+1) = material_phase(g,i,e)                    ! phaseID of grain
       c = c + 1_pInt
     case ('volume')
       crystallite_postResults(c+1) = material_volume(g,i,e)                   ! grain volume (not fraction but absolute, right?)
       c = c + 1_pInt
     case ('orientation')
       crystallite_postResults(c+1:c+4) = &
         crystallite_orientation(:,g,i,e)     ! grain orientation as quaternion
       c = c + 4_pInt
     case ('eulerangles')
       crystallite_postResults(c+1:c+3) = inDeg * & 
         math_QuaternionToEuler(crystallite_orientation(:,g,i,e)) ! grain orientation as Euler angles in degree
       c = c + 3_pInt
     case ('grainrotation')
       crystallite_postResults(c+1:c+4) = &
         math_QuaternionToAxisAngle(crystallite_rotation(1:4,g,i,e)) ! grain rotation away from initial orientation as axis-angle 
       crystallite_postResults(c+4) = inDeg * crystallite_postResults(c+4) ! angle in degree
       c = c + 4_pInt
     case ('defgrad','f')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(crystallite_partionedF(:,:,g,i,e),(/mySize/))
       c = c + mySize
     case ('fe')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(crystallite_Fe(:,:,g,i,e),(/mySize/))
       c = c + mySize
     case ('ee')
       Ee = 0.5_pReal * (math_mul33x33(transpose(crystallite_Fe(:,:,g,i,e)), crystallite_Fe(:,:,g,i,e)) - math_I3)
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(Ee(:,:),(/mySize/))
       c = c + mySize
     case ('fp')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(crystallite_Fp(:,:,g,i,e),(/mySize/))
       c = c + mySize
     case ('p','firstpiola','1stpiola')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(crystallite_P(:,:,g,i,e),(/mySize/))
       c = c + mySize
     case ('s','tstar','secondpiola','2ndpiola')
       mySize = 9_pInt
       crystallite_postResults(c+1:c+1+mySize) = reshape(math_Mandel6to33(crystallite_Tstar_v(:,g,i,e)),(/mySize/))
       c = c + mySize
   end select
 enddo
  
 crystallite_postResults(c+1) = constitutive_sizePostResults(g,i,e); c = c+1_pInt  ! size of constitutive results
 crystallite_postResults(c+1:c+constitutive_sizePostResults(g,i,e)) = &
         constitutive_postResults(crystallite_Tstar_v(:,g,i,e), crystallite_subTstar0_v(:,g,i,e), crystallite_Fe, crystallite_Fp, &
                                  crystallite_Temperature(g,i,e), crystallite_disorientation(:,:,g,i,e), dt, &
                                  crystallite_subdt(g,i,e), g, i, e)
 c = c + constitutive_sizePostResults(g,i,e)
 
 return 
 
endfunction


END MODULE
!##############################################################

