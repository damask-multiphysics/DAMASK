! SPDX-License-Identifier: AGPL-3.0-or-later
!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Utilities used by the FEM solver
!--------------------------------------------------------------------------------------------------
#include <petsc/finclude/petsc.h>
module FEM_utilities
  use PETSc
#ifndef PETSC_EXPOSES_MPI
  use MPI_f08
#endif

  use prec
  use constants
  use config
  use math
  use misc
  use IO
  use parallelization
  use discretization_mesh
  use homogenization
  use materialpoint

#ifndef PETSC_EXPOSES_MPIF90
  implicit none(type,external)
#else
  implicit none
#endif
  private

  logical,     public             :: cutBack = .false.                                              !< cut back of BVP solver in case convergence is not achieved or a material point is terminally ill

!--------------------------------------------------------------------------------------------------
! field labels information
  character(len=*), parameter, public :: &
    FIELD_MECH_label = 'mechanical'

!--------------------------------------------------------------------------------------------------
! derived types
  type, public :: tSolutionState                                                                    !< return type of solution from FEM solver variants
    logical :: converged      = .true.
    PetscInt :: iter_needed   = 0_pPETSCINT
  end type tSolutionState

  type, public :: tMechBC                                                                           !< boundary conditions data
    real(pREAL), allocatable, dimension(:) :: displacements                                         !< u_dot & u
    real(pREAL), allocatable, dimension(:) :: forces                                                !< f_dot & f
    integer,     allocatable, dimension(:) :: active                                                !< which of u_dot/u/f_dot/f is set
  end type tMechBC

  enum, bind(c); enumerator :: &                                                                    !< allowed BC types
    BC_TYPE_NONE  = int(b'0000'), &
    BC_TYPE_U_DOT = int(b'0001'), &
    BC_TYPE_U     = int(b'0010'), &
    BC_TYPE_F_DOT = int(b'0100'), &
    BC_TYPE_F     = int(b'1000')
  end enum

#if PETSC_VERSION_MINOR<23
  external :: &
    PetscSectionGetFieldComponents, &
    PetscSectionGetFieldDof, &
    PetscSectionGetFieldOffset
#endif

  public :: &
    FEM_utilities_init, &
    utilities_constitutiveResponse, &
    utilities_assembleFext, &
    utilities_assembleU, &
    utilities_projectDisplacementBC, &
    ! enum
    BC_TYPE_NONE, &
    BC_TYPE_U_DOT, &
    BC_TYPE_U , &
    BC_TYPE_F_DOT, &
    BC_TYPE_F


contains

!ToDo: use functions in variable call
!--------------------------------------------------------------------------------------------------
!> @brief Allocate all neccessary fields.
!--------------------------------------------------------------------------------------------------
subroutine FEM_utilities_init(num_mesh)

  type(tDict), pointer, intent(in) :: num_mesh                                                      !< 'mesh' dictionary from numerics file


  print'(/,1x,a)',   '<<<+-  FEM_utilities init  -+>>>'; flush(IO_STDOUT)

end subroutine FEM_utilities_init


!--------------------------------------------------------------------------------------------------
!> @brief Calculate constitutive response.
!--------------------------------------------------------------------------------------------------
subroutine utilities_constitutiveResponse(status, Delta_t,forwardData)

  integer(kind(STATUS_OK)),  intent(out)  :: status
  real(pREAL), intent(in)                 :: Delta_t                                                !< loading time
  logical,     intent(in)                 :: forwardData                                            !< age results

  integer(MPI_INTEGER_KIND) :: err_MPI


  print'(/,1x,a)', '... evaluating constitutive response ......................................'

  if (forwardData) call materialpoint_forward()
  call homogenization_mechanical_response(status,Delta_t,1,int(mesh_maxNips*mesh_nElems))           ! calculate P field
  cutBack = .false.

end subroutine utilities_constitutiveResponse


!--------------------------------------------------------------------------------------------------
!> @brief Assemble vector of external forces.
!> @details Build the end-of-step aimed external forces.
!--------------------------------------------------------------------------------------------------
subroutine utilities_assembleFext(f_aim_vec, dm, mechBC, Delta_t)

  Vec,                         intent(inout) :: f_aim_vec                                           !< (local) end-of-step aimed forces
  DM,                          intent(in)    :: dm                                                  !< (local) DM
  type(tMechBC), dimension(:), intent(in)    :: mechBC                                              !< loadstep boundary conditions data
  real(pREAL),                 intent(in)    :: Delta_t                                             !< load step total time

  PetscInt       :: component, boundary, &
                    point, dof, n_field_dof, n_field_comp, offset
  IS             :: BC_points_IS                                                                    ! IS of BC points
  PetscSection   :: section
  PetscErrorCode :: err_PETSc
  PetscInt,    pointer :: BC_points(:)                                                              ! array of IS of BC points
  real(pREAL), pointer :: f_aim(:)                                                                  ! array for f_aim_vec
  character(len=11)    :: BC_label                                                                  ! face/edge/vertex set
  integer, parameter   :: FDOT_OR_F = ior(BC_TYPE_F_DOT, BC_TYPE_F)


  ! Point loads
  call DMGetLocalSection(dm,section,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscSectionGetFieldComponents(section,0_pPETSCINT,n_field_comp,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecGetArray(f_aim_vec,f_aim,err_PETSc)
  CHKERRQ(err_PETSc)
  do boundary = 1_pPETSCINT, mesh_Nboundaries
    if (any(iand(mechBC(boundary)%active, FDOT_OR_F) > 0)) then
      BC_label = PETSC_GENERIC_LABELS(mesh_boundariesIdx(boundary))
      call DMGetStratumIS(dm,BC_label,mesh_boundariesIS(boundary),BC_points_IS,err_PETSc)
      CHKERRQ(err_PETSc)
      call ISGetIndices(BC_points_IS,BC_points,err_PETSc)
      CHKERRQ(err_PETSc)
      do point = 1_pPETSCINT, size(BC_points)
        call PetscSectionGetFieldDof(section,BC_points(point),0_pPETSCINT,n_field_dof,err_PETSc)
        CHKERRQ(err_PETSc)
        call PetscSectionGetFieldOffset(section,BC_points(point),0_pPETSCINT,offset,err_PETSc)
        CHKERRQ(err_PETSc)
        do component = 1_pPETSCINT, size(mechBC(boundary)%forces)
          if (iand(mechBC(boundary)%active(component), FDOT_OR_F) == 0) cycle
          do dof = offset+component, offset+n_field_dof, n_field_comp
            f_aim(dof) = merge(f_aim(dof) + mechBC(boundary)%forces(component)*Delta_t, &
                               mechBC(boundary)%forces(component), &
                               mechBC(boundary)%active(component) == BC_TYPE_F_DOT)
          end do
        end do
      end do
      call ISRestoreIndices(BC_points_IS,BC_points,err_PETSc)
      CHKERRQ(err_PETSc)
    end if
  end do
  call VecRestoreArray(f_aim_vec,f_aim,err_PETSc)
  CHKERRQ(err_PETSc)
  call ISDestroy(BC_points_IS,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine utilities_assembleFext


!--------------------------------------------------------------------------------------------------
!> @brief assemble vector of displacement (Dirichlet) boundary conditions
!> @details builds the end-of-step aimed displacement
!--------------------------------------------------------------------------------------------------
subroutine utilities_assembleU(u_aim_vec, dm, mechBC, Delta_t)

  Vec,                         intent(inout) :: u_aim_vec                                           !< aim displacement BC vector
  DM,                          intent(in)    :: dm                                                  !< DM (local)
  type(tMechBC), dimension(:), intent(in)    :: mechBC                                              !< BC data
  real(pREAL),                 intent(in)    :: Delta_t                                             !< load time increment

  PetscInt       :: component, boundary, &
                    point, dof, n_field_dof, n_field_comp, offset
  IS             :: BC_points_IS
  PetscSection   :: section
  PetscErrorCode :: err_PETSc
  PetscInt,    pointer :: BC_points(:)
  real(pREAL), pointer :: u_aim(:)
  character(len=11)    :: BC_label
  integer, parameter   :: UDOT_OR_U = ior(BC_TYPE_U_DOT, BC_TYPE_U)


  ! Displacement BC
  call DMGetLocalSection(dm,section,err_PETSc)
  CHKERRQ(err_PETSc)
  call PetscSectionGetFieldComponents(section,0_pPETSCINT,n_field_comp,err_PETSc)
  CHKERRQ(err_PETSc)
  call VecGetArray(u_aim_vec,u_aim,err_PETSc)
  CHKERRQ(err_PETSc)
  do boundary = 1_pPETSCINT, mesh_Nboundaries
    if (any(iand(mechBC(boundary)%active, UDOT_OR_U) > 0)) then
      BC_label = PETSC_GENERIC_LABELS(mesh_boundariesIdx(boundary))
      call DMGetStratumIS(dm,BC_label,mesh_boundariesIS(boundary),BC_points_IS,err_PETSc)
      CHKERRQ(err_PETSc)
      call ISGetIndices(BC_points_IS,BC_points,err_PETSc)
      CHKERRQ(err_PETSc)
      do point = 1_pPETSCINT, size(BC_points)
        call PetscSectionGetFieldDof(section,BC_points(point),0_pPETSCINT,n_field_dof,err_PETSc)
        CHKERRQ(err_PETSc)
        call PetscSectionGetFieldOffset(section,BC_points(point),0_pPETSCINT,offset,err_PETSc)
        CHKERRQ(err_PETSc)
        do component = 1_pPETSCINT, size(mechBC(boundary)%displacements)
          if (iand(mechBC(boundary)%active(component), UDOT_OR_U) == 0) cycle
          do dof = offset+component, offset+n_field_dof, n_field_comp
            u_aim(dof) = merge(u_aim(dof) + mechBC(boundary)%displacements(component) * Delta_t, &
                               mechBC(boundary)%displacements(component), &
                               mechBC(boundary)%active(component) == BC_TYPE_U_DOT)
          end do
        end do
      end do
      call ISRestoreIndices(BC_points_IS,BC_points,err_PETSc)
      CHKERRQ(err_PETSc)
    end if
  end do
  call VecRestoreArray(u_aim_vec,u_aim,err_PETSc)
  CHKERRQ(err_PETSc)
  call ISDestroy(BC_points_IS,err_PETSc)
  CHKERRQ(err_PETSc)

end subroutine utilities_assembleU


!--------------------------------------------------------------------------------------------------
!> @brief Add displacement boundary conditions onto current coordinates vector.
!> @details Enforce displacement BC onto the appropriate DoF by adding the end-of-step aim
!> @details displacement to current coordinates.
!--------------------------------------------------------------------------------------------------
subroutine utilities_projectDisplacementBC(x_vec, rate_vec, Delta_t)

  Vec,         intent(inout) :: x_vec
  Vec,         intent(in)    :: rate_vec
  real(preal), intent(in)    :: Delta_t

  PetscErrorCode :: err_PETSc


  call VecAXPY(x_vec, Delta_t, rate_vec, err_PETSc)                                                 ! x = x + rate * dt
  CHKERRQ(err_PETSc)

end subroutine utilities_projectDisplacementBC

end module FEM_utilities
