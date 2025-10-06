!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Driver controlling inner and outer load case looping of the FEM solver
!> @details doing cutbacking, reporting statistics, writing
!> results
!--------------------------------------------------------------------------------------------------
program DAMASK_mesh
#include <petsc/finclude/petscsys.h>
  use PetscDM
  use prec
  use CLI
  use parallelization
  use materialpoint
  use discretization_mesh
  use FEM_Utilities
  use mesh_mechanical_FEM

  implicit none(type,external)

  type :: tLoadCase
    real(pREAL)   :: t                   = 0.0_pREAL                                                !< length of increment
    integer       :: N                   = 0, &                                                     !< number of increments
                     f_out               = 1                                                        !< frequency of result writes
    logical       :: estimate_rate       = .true.                                                   !< follow trajectory of former loadcase
    type(tMechBC),  allocatable, dimension(:) :: mechBC
  end type tLoadCase


!--------------------------------------------------------------------------------------------------
! loop variables, convergence etc.
  integer, parameter :: &
    subStepFactor = 2                                                                               !< for each substep, divide the last time increment by 2.0
  real(pREAL) :: &
    t   = 0.0_pREAL, &                                                                              !< elapsed time
    t_0 = 0.0_pREAL, &                                                                              !< begin of interval
    Delta_t = 0.0_pREAL, &                                                                          !< current time interval
    Delta_t_prev = 0.0_pREAL                                                                        !< previous time interval
  logical :: &
    guess, &                                                                                        !< guess along former trajectory
    stagIterate
  integer :: &
    l, &
    m, &
    errorID, &
    cutBackLevel = 0, &                                                                             !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
    stepFraction = 0, &                                                                             !< fraction of current time interval
    inc, &                                                                                          !< current increment in current load case
    totalIncsCounter = 0, &                                                                         !< total # of increments
    statUnit = 0, &                                                                                 !< file unit for statistics output
    stagIter, &
    component
  type(tDict), pointer :: &
    num_solver, &
    num_mesh, &
    load, &
    load_step, &
    step_bc, &
    mech_BC, &
    step_discretization
  type(tList), pointer :: &
    load_steps, &
    mech_u, &
    step_mech
  character(len=pSTRLEN) :: &
    incInfo, &
    loadcase_string
  integer :: &
    stagItMax, &                                                                                    !< max number of field level staggered iterations
    maxCutBack, &                                                                                   !< max number of cutbacks
    skipBCTagDigits, &                                                                              !< number of characters to skip to print BC tag (T descriptor)
    BCTag                                                                                           !< tag value read from YAML load file
  integer, allocatable, dimension(:) :: &
    knownBCTags                                                                                     !< array of read BC tags (check for duplicates)


  type(tLoadCase), allocatable, dimension(:) :: loadCases                                           !< array of all load cases
  type(tSolutionState), allocatable, dimension(:) :: solres
  PetscInt :: boundary, dimPlex
  PetscErrorCode :: err_PETSc
  external :: &
    quit
  character(len=:), allocatable :: &
    fileContent, fname, tagPrintFormat
  character(len=6) :: BC_elem


!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
  call materialpoint_initAll()
  print'(/,1x,a)', '<<<+-  DAMASK_mesh init  -+>>>'; flush(IO_STDOUT)

!---------------------------------------------------------------------
! reading field information from numerics file and do sanity checks
  num_solver => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_mesh   => num_solver%get_dict('mesh',defaultVal=emptyDict)
  stagItMax  = num_mesh%get_asInt('N_staggered_iter_max',defaultVal=10)
  maxCutBack = num_mesh%get_asInt('N_cutback_max',defaultVal=3)

  if (stagItMax < 0)  call IO_error(301,ext_msg='N_staggered_iter_max')
  if (maxCutBack < 0) call IO_error(301,ext_msg='N_cutback_max')

! reading basic information from load case file and allocate data structure containing load cases
  call DMGetDimension(geomMesh,dimPlex,err_PETSc)                                                   !< dimension of mesh (2D or 3D)
  CHKERRA(err_PETSc)
  allocate(solres(1))

  if (worldrank == 0) then
    fileContent = IO_read(CLI_loadFile)
    fname = CLI_loadFile
    if (scan(fname,'/') /= 0) fname = fname(scan(fname,'/',.true.)+1:)
    call result_openJobFile(parallel=.false.)
    call result_addSetupFile(fileContent,fname,'load case definition (mesh solver)')
    call result_closeJobFile()
  end if

  call parallelization_bcast_str(fileContent)
  load => YAML_str_asDict(fileContent)
  load_steps => load%get_list('loadstep')

  allocate(loadCases(size(load_steps)))
  do l = 1, size(load_steps)
    load_step => load_steps%get_dict(l)
    step_bc   => load_step%get_dict('boundary_conditions')
    step_mech => step_bc%get_list('mechanical')
    allocate(loadCases(l)%mechBC(mesh_Nboundaries))
    loadCases(l)%mechBC(:)%nComponents = dimPlex                                                   !< X, Y (, Z) displacements
    do boundary = 1_pPETSCINT, mesh_Nboundaries
      allocate(loadCases(l)%mechBC(boundary)%dot_u(dimPlex),  source = 0.0_pREAL)
      allocate(loadCases(l)%mechBC(boundary)%active(dimPlex), source = .false.)
    end do

    allocate(knownBCTags(size(step_mech)), source = -1)
    do m = 1, size(step_mech)
      mech_BC => step_mech%get_dict(m)
      BCTag = mech_BC%get_asInt('tag')
      boundary = findloc(mesh_boundariesIS, BCtag, dim = 1)
      if (boundary == 0) then                                                                      !< tag not defined in mesh file
        call IO_error(error_ID = 837, ext_msg = 'BC tag '//mech_BC%get_asStr('tag')// &
                                                ' refers to inexistent node/edge/face mesh group')
      else if (findloc(knownBCTags, BCTag, dim = 1) /= 0) then                                     !< duplicated tag
        call IO_error(error_ID = 837, ext_msg = 'BC redefinition not allowed (tag '// &
                                                mech_BC%get_asStr('tag')//')')
      end if
      knownBCTags(m) = BCTag
      mech_u => mech_BC%get_list('dot_u')
      do component = 1, dimPlex
        if (mech_u%get_asStr(component) /= 'x') then
          loadCases(l)%mechBC(boundary)%active(component) = .true.
          loadCases(l)%mechBC(boundary)%dot_u(component)  = mech_u%get_asReal(component)
        end if
      end do
    end do
    deallocate(knownBCTags)

    step_discretization => load_step%get_dict('discretization')
    loadCases(l)%t = step_discretization%get_asReal('t')
    loadCases(l)%N = step_discretization%get_asInt  ('N')

    if (load_step%get_asStr('f_out',defaultVal='n/a') == 'none') then
      loadCases(l)%f_out = huge(0)
    else
      loadCases(l)%f_out = load_step%get_asInt('f_out', defaultVal=1)
    end if
    loadCases(l)%estimate_rate = (load_step%get_asBool('estimate_rate',defaultVal=.true.) .and. l>1)
  end do

!--------------------------------------------------------------------------------------------------
! consistency checks and output of load case
  errorID = 0
  skipBCTagDigits = 4+6+1+floor(log10(real(maxval(mesh_boundariesIS))))+1+2                        !< Indentation(4)+BC_elem(6)+blank(1)+descriptor T(1)+NumDigits(floor(..)+1)+2(start printing after 1 blank)
  allocate(character(len=skipBCTagDigits) :: tagPrintFormat)
  write(tagPrintFormat,'(i0)') skipBCTagDigits
  checkLoadcases: do l = 1, size(load_steps)
    write (loadcase_string, '(i0)' ) l
    print'(/,1x,a,1x,i0)', 'load case:', l
    if (.not. loadCases(l)%estimate_rate) &
      print'(2x,a)', 'drop guessing along trajectory'
    print'(2x,a)', 'Field '//trim(FIELD_MECH_label)

    do boundary = 1_pPETSCINT, mesh_Nboundaries
      BC_elem = merge('Vertex', merge('Edge  ', 'Face  ', dimplex == 2_pPETSCINT), &
                      mesh_boundariesIdx(boundary) == mesh_BCTypeVertex)
      do component = 1_pPETSCINT, dimPlex
        if (loadCases(l)%mechBC(boundary)%active(component)) &
          print'(4x,a,1x,i0,T'//tagPrintFormat//',a,1x,i1,1x,a,1x,f12.7)', &
          BC_elem, mesh_boundariesIS(boundary), &
          'Component', component, &
          'Value', loadCases(l)%mechBC(boundary)%dot_u(component)
      end do
    end do

    print'(2x,a,T21,g0.6)', 'time:', loadCases(l)%t
    if (loadCases(l)%N < 1) errorID = 835                                                          ! non-positive incs count
    print'(2x,a,T21,i0)',   'increments:', loadCases(l)%N
    if (loadCases(l)%f_out < 1) errorID = 836                                                      ! non-positive result frequency
    print'(2x,a,T21,i0)',   'output frequency:', loadCases(l)%f_out
    if (errorID > 0) call IO_error(error_ID = errorID, ext_msg = loadcase_string)                  ! exit with error message
  end do checkLoadcases

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call FEM_Utilities_init(num_mesh)
  call FEM_mechanical_init(loadCases(1)%mechBC,num_mesh)
  call config_numerics_deallocate()

  if (worldrank == 0) then
    open(newunit=statUnit,file=trim(CLI_jobName)//'.sta',form='FORMATTED',status='REPLACE')
    write(statUnit,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                  ! statistics file
  end if

  print'(/,1x,a)', '... saving initial configuration ..........................................'
  flush(IO_STDOUT)
  call materialpoint_result(0,0.0_pREAL)

  loadCaseLooping: do l = 1, size(load_steps)
    t_0 = t                                                                                        ! load case start time
    guess = loadCases(l)%estimate_rate                                                             ! change of load case? homogeneous guess for the first inc

    incLooping: do inc = 1, loadCases(l)%N
      totalIncsCounter = totalIncsCounter + 1

!--------------------------------------------------------------------------------------------------
! forwarding time
      Delta_t_prev = Delta_t                                                                        ! last timeinc that brought former inc to an end
      Delta_t = loadCases(l)%t/real(loadCases(l)%N,pREAL)
      Delta_t = Delta_t * real(subStepFactor,pREAL)**real(-cutBackLevel,pREAL)                      ! depending on cut back level, decrease time step
      stepFraction = 0                                                                              ! fraction scaled by stepFactor**cutLevel

      subStepLooping: do while (stepFraction < subStepFactor**cutBackLevel)
        t = t + Delta_t                                                                             ! forward target time
        stepFraction = stepFraction + 1                                                             ! count step

!--------------------------------------------------------------------------------------------------
! report begin of new step
        print'(/,1x,a)', '###########################################################################'
        print'(1x,a,es12.5,6(a,i0))',&
                'Time', t, &
                's: Increment ', inc, '/', loadCases(l)%N,&
                '-', stepFraction, '/', subStepFactor**cutBackLevel,&
                ' of load case ', l,'/', size(load_steps)
        write(incInfo,'(4(a,i0))') &
               'Increment ',totalIncsCounter,'/',sum(loadCases%N),&
               '-',stepFraction, '/', subStepFactor**cutBackLevel
        flush(IO_STDOUT)

        call FEM_mechanical_forward(guess,Delta_t,Delta_t_prev,loadCases(l)%mechBC)

!--------------------------------------------------------------------------------------------------
! solve fields
        stagIter = 0
        stagIterate = .true.
        do while (stagIterate)
          solres(1) = FEM_mechanical_solution(incInfo,Delta_t,Delta_t_prev,loadCases(l)%mechBC)
          if (.not. solres(1)%converged) exit

          stagIter = stagIter + 1
          stagIterate =            stagIter < stagItMax &
                       .and.       all(solres(:)%converged) &
                       .and. .not. all(solres(:)%stagConverged)                                     ! stationary with respect to staggered iteration
        end do

! check solution
        cutBack = .False.
        if (.not. all(solres(:)%converged .and. solres(:)%stagConverged)) then                      ! no solution found
          if (cutBackLevel < maxCutBack) then                                                       ! do cut back
            cutBack = .True.
            stepFraction = (stepFraction - 1) * subStepFactor                                       ! adjust to new denominator
            cutBackLevel = cutBackLevel + 1
            t    = t - Delta_t                                                                      ! rewind time
            Delta_t = Delta_t/2.0_pREAL
            print'(/,1x,a)', 'cutting back'
          else                                                                                      ! default behavior, exit if spectral solver does not converge
            if (worldrank == 0) close(statUnit)
            call IO_error(950)
          end if
        else
          guess = .true.                                                                            ! start guessing after first converged (sub)inc
          Delta_t_prev = Delta_t
        end if
        if (.not. cutBack .and. worldrank == 0) then
          write(statUnit,*) totalIncsCounter, t, cutBackLevel, &
                            solres%converged, solres%iterationsNeeded                               ! write statistics about accepted solution
          flush(statUnit)
        end if
      end do subStepLooping

      cutBackLevel = max(0, cutBackLevel - 1)                                                       ! try half number of subincs next inc

      if (all(solres(:)%converged)) then
        print'(/,1x,a,1x,i0,1x,a)', 'increment', totalIncsCounter, 'converged'
      else
        print'(/,1x,a,1x,i0,1x,a)', 'increment', totalIncsCounter, 'NOT converged'
      end if; flush(IO_STDOUT)

      if (mod(inc,loadCases(l)%f_out) == 0) then                            ! at output frequency
        print'(/,1x,a)', '... saving results ........................................................'
        call FEM_mechanical_updateCoords()
        call materialpoint_result(totalIncsCounter,t)
      end if


    end do incLooping

  end do loadCaseLooping


!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  print'(/,1x,a)', '###########################################################################'
  if (worldrank == 0) close(statUnit)

  call quit(0)                                                                                      ! no complains ;)


end program DAMASK_mesh
