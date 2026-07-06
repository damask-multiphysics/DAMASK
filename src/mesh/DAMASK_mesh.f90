! SPDX-License-Identifier: AGPL-3.0-or-later
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
  use IO
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
    sub_step_factor = 2                                                                             !< for each substep, divide the last time increment by 2.0
  real(pREAL) :: &
    t   = 0.0_pREAL, &                                                                              !< elapsed time
    t_0 = 0.0_pREAL, &                                                                              !< begin of interval
    Delta_t = 0.0_pREAL, &                                                                          !< current time interval
    Delta_t_prev = 0.0_pREAL                                                                        !< previous time interval
  logical :: &
    guess
  integer :: &
    l, &
    cut_back_level = 0, &                                                                           !< cut back level \f$ t = \frac{t_{inc}}{2^l} \f$
    step_fraction = 0, &                                                                            !< fraction of current time interval
    inc, &                                                                                          !< current increment in current load case
    n_total_inc = 0, &                                                                              !< total # of increments
    unit_stat = 0                                                                                   !< file unit for statistics output
  type(tDict), pointer :: &
    num_solver => NULL(), &                                                                         !< dictionaries: solver
    num_mesh => NULL(), &                                                                           !<               mechanical
    load => NULL()                                                                                  !<               full load file
  character(len=pSTRLEN) :: &
    inc_info                                                                                        !< report new step (iterations, cut backs, etc)
  integer :: &
    max_cutback                                                                                     !< max number of cutbacks
  type(tLoadCase), allocatable, dimension(:) :: &
    load_cases                                                                                      !< array of all load cases
  type(tSolutionState) :: &
    sol_state                                                                                       !< solution status (converged, iters needed)
  character(len=:), allocatable :: &
    f_content, f_name                                                                               !< load file full content/filename


!--------------------------------------------------------------------------------------------------
! init DAMASK (all modules)
  call materialpoint_initAll()
  print'(/,1x,a)', '<<<+-  DAMASK_mesh init  -+>>>'; flush(IO_STDOUT)

!--------------------------------------------------------------------------------------------------
! reading field information from numerics file and do sanity checks
  num_solver => config_numerics%get_dict('solver',defaultVal=emptyDict)
  num_mesh   => num_solver%get_dict('mesh',defaultVal=emptyDict)
  max_cutback   = num_mesh%get_asInt('N_cutback_max',defaultVal=3)

  if (max_cutback < 0) call IO_error(301,ext_msg='N_cutback_max')

!--------------------------------------------------------------------------------------------------
! parse load life, do all necessary checks for correct BC, and output load case information
  if (worldrank == 0) then
    f_content = IO_read(CLI_loadFile)
    f_name = CLI_loadFile
    if (scan(f_name,'/') /= 0) f_name = f_name(scan(f_name,'/',.true.)+1:)
    call result_openJobFile(parallel=.false.)
    call result_addSetupFile(f_content,f_name,'load case definition (mesh solver)')
    call result_closeJobFile()
  end if

  call parallelization_bcast_str(f_content)
  load => YAML_str_asDict(f_content)
  load_cases = parse_and_print_load_cases(load%get_list('loadstep'))

!--------------------------------------------------------------------------------------------------
! doing initialization depending on active solvers
  call FEM_Utilities_init(num_mesh)
  call FEM_mechanical_init(load_cases(1)%mechBC,num_mesh)
  call config_numerics_deallocate()

  if (worldrank == 0) then
    open(newunit=unit_stat,file=trim(CLI_jobName)//'.sta',form='FORMATTED',status='REPLACE')
    write(unit_stat,'(a)') 'Increment Time CutbackLevel Converged IterationsNeeded'                 ! statistics file
  end if

  print'(/,1x,a)', '... saving initial configuration ..........................................'
  flush(IO_STDOUT)
  call materialpoint_result(0,0.0_pREAL)

  loadCaseLooping: do l = 1, size(load_cases)
    t_0 = t                                                                                         ! load case start time
    guess = load_cases(l)%estimate_rate                                                             ! change of load case? homogeneous guess for the first inc
    call FEM_mechanical_assembleFext(load_cases(l)%mechBC, load_cases(l)%t)                         ! assemble external loads vector
    call FEM_mechanical_assembleU(load_cases(l)%mechBC, load_cases(l)%t)                            ! assemble vector of displacements (Dirichlet) BC

    incLooping: do inc = 1, load_cases(l)%N
      n_total_inc = n_total_inc + 1

!--------------------------------------------------------------------------------------------------
! forwarding time
      Delta_t_prev = Delta_t                                                                        ! last timeinc that brought former inc to an end
      Delta_t = load_cases(l)%t/real(load_cases(l)%N,pREAL)
      Delta_t = Delta_t * real(sub_step_factor,pREAL)**real(-cut_back_level,pREAL)                  ! depending on cut back level, decrease time step
      step_fraction = 0                                                                             ! fraction scaled by stepFactor**cutLevel

      subStepLooping: do while (step_fraction < sub_step_factor**cut_back_level)
        t = t + Delta_t                                                                             ! forward target time
        step_fraction = step_fraction + 1                                                           ! count step

!--------------------------------------------------------------------------------------------------
! report begin of new step
        print'(/,1x,a)', '###########################################################################'
        print'(1x,a,es12.5,6(a,i0))',&
                'Time', t, &
                's: Increment ', inc, '/', load_cases(l)%N,&
                '-', step_fraction, '/', sub_step_factor**cut_back_level,&
                ' of load case ', l,'/', size(load_cases)
        write(inc_info,'(4(a,i0))') &
               'Increment ',n_total_inc,'/',sum(load_cases%N),&
               '-',step_fraction, '/', sub_step_factor**cut_back_level
        flush(IO_STDOUT)

        call FEM_mechanical_forward(guess,Delta_t,Delta_t_prev)

!--------------------------------------------------------------------------------------------------
! solve fields
        sol_state = FEM_mechanical_solution(inc_info,Delta_t,Delta_t_prev,load_cases(l)%mechBC)
        cutBack = .False.

        if (.not. sol_state%converged) then                                                         ! no solution found
          if (cut_back_level < max_cutback) then                                                    ! do cut back
            cutBack = .True.
            step_fraction = (step_fraction - 1) * sub_step_factor                                   ! adjust to new denominator
            cut_back_level = cut_back_level + 1
            t = t - Delta_t                                                                         ! rewind time
            Delta_t = Delta_t/2.0_pREAL
            print'(/,1x,a)', 'cutting back'
          else                                                                                      ! default behavior, exit if spectral solver does not converge
            if (worldrank == 0) close(unit_stat)
            call IO_error(950)
          end if
        else
          guess = .true.                                                                            ! start guessing after first converged (sub)inc
          Delta_t_prev = Delta_t
        end if
        if (.not. cutBack .and. worldrank == 0) then
          write(unit_stat,*) n_total_inc, t, cut_back_level, &
                            sol_state%converged, sol_state%iter_needed                              ! write statistics about accepted solution
          flush(unit_stat)
        end if
      end do subStepLooping

      cut_back_level = max(0, cut_back_level - 1)                                                   ! try half number of subincs next inc

      if (mod(inc,load_cases(l)%f_out) == 0) then                                                   ! at output frequency
        print'(/,1x,a)', '... saving results ........................................................'
        call FEM_mechanical_updateCoords()
        call materialpoint_result(n_total_inc,t)
      end if
    end do incLooping
  end do loadCaseLooping

!--------------------------------------------------------------------------------------------------
! report summary of whole calculation
  print'(/,1x,a)', '###########################################################################'
  if (worldrank == 0) close(unit_stat)

  call quit(0)                                                                                      ! no complains ;)


contains

!--------------------------------------------------------------------------------------------------
!> @brief Parse load cases data from YAML file, perform all necessary checks, and print it.
!--------------------------------------------------------------------------------------------------
function parse_and_print_load_cases(load_steps) result(load_cases)

#ifndef __GFORTRAN__
#include <petsc/finclude/petscsys.h>
  use PetscDM

  use prec
  use types
  use IO
  use discretization_mesh
  use FEM_Utilities
  use mesh_mechanical_FEM

  import, only: tLoadCase
#endif

  type(tList), target, intent(in) :: load_steps                                                     !< full YAML file dictionary
  type(tLoadCase), allocatable, dimension(:) :: load_cases                                          !< array of all load cases

  integer :: &
    boundary, component, l, m, &                                                                    ! loop variables
    BC_tag                                                                                          ! tag used for BC
  PetscInt :: &
    dimPlex                                                                                         ! mesh dimension (2D or 3D)
  real(pREAL) :: &
    BC_comp_value                                                                                   ! BC x/y/z component value from YAML load file
  logical, dimension(mesh_nBoundaries) :: &
    known_BC, &                                                                                     ! already read BC entries (repetition per load-step forbidden)
    use_label                                                                                       ! BC uses label or tag
  character(len=pSTRLEN) :: &
    BC_label                                                                                        ! mesh label or generic label
  character(len=3) :: &
    BC_unit                                                                                         ! unit for printing: one of [m, m/s, N, N/s]
  PetscErrorCode :: &
    err_PETSc
  type(tDict), pointer :: &
    load_step           => NULL(), &                                                                ! dictionaries: loadstep entry
    boundary_conditions => NULL(), &                                                                !               boundary conditions entry
    BC_mechanical       => NULL(), &                                                                !               single BC entry (mechanics solver)
    step_discretization => NULL()                                                                   !               discretization entry
  type(tList), pointer :: &
    BCs_mechanical      => NULL(), &                                                                !               mechanical entry
    BC_load_comps       => NULL()                                                                   ! u/f[_dot] components list


  call DMGetDimension(geomMesh, dimPlex, err_PETSc)
  CHKERRA(err_PETSc)

!--------------------------------------------------------------------------------------------------
! allocate loadstep-related arrays
  allocate(load_cases(size(load_steps)))
  do l = 1, size(load_steps)
    load_step => load_steps%get_dict(l)
    boundary_conditions   => load_step%get_dict('boundary_conditions')
    BCs_mechanical => boundary_conditions%get_list('mechanical')
    allocate(load_cases(l)%mechBC(mesh_nBoundaries))
    do boundary = 1, int(mesh_nBoundaries)
      allocate(load_cases(l)%mechBC(boundary)%displacements(dimPlex), source = 0.0_pREAL)
      allocate(load_cases(l)%mechBC(boundary)%forces(dimPlex), source = 0.0_pREAL)
      allocate(load_cases(l)%mechBC(boundary)%active(dimPlex), source = BC_TYPE_NONE)
    end do

!--------------------------------------------------------------------------------------------------
! check valid tags/labels
    use_label = .false.                                                                             ! reset per load step
    known_BC = .false.
    do m = 1, size(BCs_mechanical)
      BC_mechanical => BCs_mechanical%get_dict(m)
      if (BC_mechanical%contains('label') .and. BC_mechanical%contains('tag')) then
        call IO_error(812_pI16, '"label" and "tag" are given for boundary condition', m, emph=[2])
      elseif (BC_mechanical%contains('label')) then
        BC_label = BC_mechanical%get_asStr('label')
        boundary = findloc(mesh_BCLabels, BC_label, dim = 1)
        if (boundary == 0) &                                                                        ! label not defined in mesh file
          call IO_error(812_pI16, 'label', trim(BC_label), 'not defined', emph = [2])
        if (known_BC(boundary)) &                                                                   ! duplicated label/tag
          call IO_error(812_pI16, 'duplicated entries: label', trim(BC_label), 'and tag', &
                        mesh_boundariesIS(boundary), emph = [2,4])
        known_BC(boundary) = .true.
        use_label(boundary) = .true.
      else if (BC_mechanical%contains('tag')) then
        BC_tag = BC_mechanical%get_asInt('tag')
        boundary = findloc(mesh_boundariesIS, BC_tag, dim = 1)
        if (boundary == 0) &                                                                        ! tag not defined in mesh file
          call IO_error(812_pI16, 'tag', BC_tag, 'not defined', emph = [2])
        if (known_BC(boundary)) &                                                                   ! duplicated tag/label
          call IO_error(812_pI16, 'duplicated entries: tag', BC_tag, 'and label', &
                        trim(mesh_BCLabels(boundary)), emph = [2, 4])
        known_BC(boundary) = .true.
      else
        call IO_error(812_pI16, 'neither "label" nor "tag" given for boundary condition', m, emph=[2])
      end if
!--------------------------------------------------------------------------------------------------
! read BC data
      associate (BC => load_cases(l)%mechBC(boundary))
        if (BC_mechanical%contains('u_dot')) then
          BC_load_comps => BC_mechanical%get_list('u_dot')
        else if (BC_mechanical%contains('dot_u')) then
          BC_load_comps => BC_mechanical%get_list('dot_u')
        end if
        if (associated(BC_load_comps)) then
          do component = 1, int(dimPlex)
            if (BC_load_comps%get_asStr(component) /= 'x') then
              BC%displacements(component) = BC_load_comps%get_asReal(component)
              BC%active(component) = ior(BC%active(component), BC_TYPE_U_DOT)
            end if
          end do
          nullify(BC_load_comps)
        end if
        if (BC_mechanical%contains('u')) then
          BC_load_comps => BC_mechanical%get_list('u')
          do component = 1, int(dimPlex)
            if (BC_load_comps%get_asStr(component) /= 'x') then
              BC%displacements(component) = BC_load_comps%get_asReal(component)
              BC%active(component) = ior(BC%active(component), BC_TYPE_U)
            end if
          end do
          nullify(BC_load_comps)
        end if
        if (BC_mechanical%contains('f_dot')) then
          BC_load_comps => BC_mechanical%get_list('f_dot')
        else if (BC_mechanical%contains('dot_f')) then
          BC_load_comps => BC_mechanical%get_list('dot_f')
        end if
        if (associated(BC_load_comps)) then
          do component = 1, int(dimPlex)
            if (BC_load_comps%get_asStr(component) /= 'x') then
              BC%forces(component) = BC_load_comps%get_asReal(component)
              BC%active(component) = ior(BC%active(component), BC_TYPE_F_DOT)
            end if
          end do
          nullify(BC_load_comps)
        end if
        if (BC_mechanical%contains('f')) then
          BC_load_comps => BC_mechanical%get_list('f')
          do component = 1, int(dimPlex)
            if (BC_load_comps%get_asStr(component) /= 'x') then
              BC%forces(component) = BC_load_comps%get_asReal(component)
              BC%active(component) = ior(BC%active(component), BC_TYPE_F)
            end if
          end do
          nullify(BC_load_comps)
        end if
!--------------------------------------------------------------------------------------------------
! check valid BC definition
        do component = 1, int(dimPlex)
          if (popcnt(BC%active(component)) > 1) then
            if (use_label(boundary)) then
              call IO_error(812_pI16, 'more than one condition specified for component', &
                            component, IO_EOL, 'in loadcase', l, '/ label', trim(BC_label), emph = [2,5,7])
            else
              call IO_error(812_pI16, 'more than one condition specified for component', &
                            component, IO_EOL, 'in loadcase', l, '/ tag', BC_tag, emph = [2,5,7])
            end if
          end if
        end do
      end associate
    end do

!--------------------------------------------------------------------------------------------------
! store discretization, time and frequency; check values
    step_discretization => load_step%get_dict('discretization')
    load_cases(l)%t = step_discretization%get_asReal('t')
    if (load_cases(l)%t < 0.0_pREAL) &
      call IO_error(301_pI16, 'loadcase', l, 'has non-positive time step length', 't', emph = [2,4])

    load_cases(l)%N = step_discretization%get_asInt('N')
    if (load_cases(l)%N < 1) &
      call IO_error(301_pI16, 'loadcase', l, 'has non-positive number of steps', 'N', emph = [2,4])

    if (load_step%get_asStr('f_out',defaultVal='n/a') == 'none') then
      load_cases(l)%f_out = huge(0)
    else
      load_cases(l)%f_out = load_step%get_asInt('f_out', defaultVal=1)
    end if
    if (load_cases(l)%f_out < 1) &
      call IO_error(301_pI16, 'loadcase', l, 'has non-positive output frequency', 'f_out', emph = [2,4])

    load_cases(l)%estimate_rate = (load_step%get_asBool('estimate_rate',defaultVal=.true.) .and. l>1)

!--------------------------------------------------------------------------------------------------
! output of load case information
    print'(/,1x,a,1x,i0)', 'load case:', l
    if (.not. load_cases(l)%estimate_rate) print'(2x,a)', 'drop guessing along trajectory'
    print'(2x,a)', 'Field '//trim(FIELD_MECH_label)

    do boundary = 1, int(mesh_nBoundaries)
      associate (BC => load_cases(l)%mechBC(boundary))
        if (all(BC%active == BC_TYPE_NONE)) cycle
        if (use_label(boundary)) then
          BC_label = mesh_BCLabels(boundary)
        else
          m = mesh_boundariesIdx(boundary)
          if (dimPlex == 2_pPETSCINT .and. m < size(PETSC_GENERIC_LABELS)) m = m + 1                ! adjust for 2D (cells -> faces)
          BC_label = PETSC_GENERIC_LABELS(m)
        end if

        print'(3x,a,1x,a,i0,a)', &
          trim(BC_label), '(', mesh_boundariesIS(boundary), ')'
        do component = 1_pPETSCINT, dimPlex
          if (BC%active(component) == BC_TYPE_NONE) cycle
          select case (BC%active(component))
            case (BC_TYPE_U_DOT : BC_TYPE_U)
              BC_comp_value = BC%displacements(component)
              BC_unit = merge('m/s', 'm  ', BC%active(component) == BC_TYPE_U_DOT)
            case (BC_TYPE_F_DOT : BC_TYPE_F)
              BC_comp_value = BC%forces(component)
              BC_unit = merge('N/s', 'N  ', BC%active(component) == BC_TYPE_F_DOT)
            end select
          print'(5x,a,1x,i1,a,1x,en12.3e2,2x,a)', &
            'Component', component, ':', BC_comp_value, BC_unit
        end do
      end associate
    end do
    print'(2x,a,T19,en12.3e2,2x,a)', 'time:',             load_cases(l)%t, 's'
    print'(2x,a,T22,i0)',            'increments:',       load_cases(l)%N
    print'(2x,a,T22,i0)',            'output frequency:', load_cases(l)%f_out
  end do

end function parse_and_print_load_cases

end program DAMASK_mesh
