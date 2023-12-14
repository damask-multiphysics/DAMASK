!--------------------------------------------------------------------------------------------------
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Luc Hantcherli, Max-Planck-Institut für Eisenforschung GmbH
!> @author W.A. Counts
!> @author Denny Tjahjanto, Max-Planck-Institut für Eisenforschung GmbH
!> @author Christoph Kords, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Interfaces DAMASK with MSC.Marc
!--------------------------------------------------------------------------------------------------
#define QUOTE(x) #x
#define PASTE(x,y) x ## y

#ifdef Marc4DAMASK
#define MARC4DAMASK Marc4DAMASK
#endif

#include "../prec.f90"
#include "../constants.f90"
#include "../parallelization.f90"
#include "../misc.f90"
#include "../IO.f90"
#include "../types.f90"
#include "../YAML_parse.f90"
#include "../HDF5_utilities.f90"

module DAMASK_interface
  use, intrinsic :: ISO_fortran_env, only: &
    compiler_version, &
    compiler_options
  use ifport, only: &
    CHDIR

  use prec
  use IO

  implicit none(type,external)
  private

  logical,          protected, public :: symmetricSolver
  character(len=*), parameter, public :: INPUTFILEEXTENSION = '.dat'

  public :: &
    DAMASK_interface_init, &
    getSolverJobName

contains

!--------------------------------------------------------------------------------------------------
!> @brief reports and sets working directory
!--------------------------------------------------------------------------------------------------
subroutine DAMASK_interface_init

  integer, dimension(8)   :: dateAndTime
  integer                 :: ierr
  character(len=pPathLen) :: wd

  external                :: quit

  print'(/,1x,a)', '<<<+-  DAMASK_Marc init -+>>>'

  print*, 'Roters et al., Computational Materials Science 158:420–478, 2019'
  print*, 'https://doi.org/10.1016/j.commatsci.2018.04.030'

  print'(/,a)', ' Version: '//DAMASKVERSION

  print'(/,a)', ' Compiled with: '//compiler_version()
  print'(a)',   ' Compiler options: '//compiler_options()

  print'(/,a)', ' Compiled on: '//__DATE__//' at '//__TIME__

  call date_and_time(values = dateAndTime)
  print'(/,a,2(i2.2,a),i4.4)', ' Date: ',dateAndTime(3),'/',dateAndTime(2),'/', dateAndTime(1)
  print'(a,2(i2.2,a),i2.2)',   ' Time: ',dateAndTime(5),':', dateAndTime(6),':', dateAndTime(7)

  inquire(5, name=wd)
  wd = wd(1:scan(wd,'/',back=.true.))
  ierr = CHDIR(wd)
  if (ierr /= 0) then
    print*, 'working directory "'//trim(wd)//'" does not exist'
    call quit(1)
  end if
  symmetricSolver = solverIsSymmetric()

end subroutine DAMASK_interface_init


!--------------------------------------------------------------------------------------------------
!> @brief solver job name (no extension) as combination of geometry and load case name
!--------------------------------------------------------------------------------------------------
function getSolverJobName()

  character(len=:), allocatable :: getSolverJobName
  character(1024)               :: inputName
  character(len=*), parameter   :: pathSep = achar(47)//achar(92)                                   ! forward and backward slash
  integer :: extPos

  inquire(5, name=inputName)                                                                        ! determine inputfile
  extPos = len_trim(inputName)-4
  getSolverJobName=inputName(scan(inputName,pathSep,back=.true.)+1:extPos)

end function getSolverJobName


!--------------------------------------------------------------------------------------------------
!> @brief determines whether a symmetric solver is used
!--------------------------------------------------------------------------------------------------
logical function solverIsSymmetric()

  character(len=pSTRLEN) :: line
  integer :: myStat,fileUnit,s,e

  open(newunit=fileUnit, file=getSolverJobName()//INPUTFILEEXTENSION, &
       status='old', position='rewind', action='read',iostat=myStat)
  do
    read (fileUnit,'(A)',END=100) line
    if (index(trim(IO_lc(line)),'solver') == 1) then
      read (fileUnit,'(A)',END=100) line                                                            ! next line
        s =     verify(line,      ' ')                                                              ! start of first chunk
        s = s + verify(line(s+1:),' ')                                                              ! start of second chunk
        e = s + scan  (line(s+1:),' ')                                                              ! end of second chunk
      solverIsSymmetric = line(s:e) /= '1'
    end if
  end do
100 close(fileUnit)

end function solverIsSymmetric

end module DAMASK_interface

#include "../result.f90"
#include "../config.f90"
#include "../LAPACK_interface.f90"
#include "../math.f90"
#include "../rotations.f90"
#include "../polynomials.f90"
#include "../tables.f90"
#include "../crystal.f90"
#include "element.f90"
#include "../geometry_plastic_nonlocal.f90"
#include "../discretization.f90"
#include "discretization_Marc.f90"
#include "../material.f90"
#include "../phase.f90"
#include "../phase_mechanical.f90"
#include "../phase_mechanical_elastic.f90"
#include "../phase_mechanical_plastic.f90"
#include "../phase_mechanical_plastic_none.f90"
#include "../phase_mechanical_plastic_isotropic.f90"
#include "../phase_mechanical_plastic_phenopowerlaw.f90"
#include "../phase_mechanical_plastic_kinehardening.f90"
#include "../phase_mechanical_plastic_dislotwin.f90"
#include "../phase_mechanical_plastic_dislotungsten.f90"
#include "../phase_mechanical_plastic_nonlocal.f90"
#include "../phase_mechanical_eigen.f90"
#include "../phase_mechanical_eigen_thermalexpansion.f90"
#include "../phase_thermal.f90"
#include "../phase_thermal_source_dissipation.f90"
#include "../phase_thermal_source_externalheat.f90"
#include "../phase_damage.f90"
#include "../phase_damage_isobrittle.f90"
#include "../phase_damage_anisobrittle.f90"
#include "../homogenization.f90"
#include "../homogenization_mechanical.f90"
#include "../homogenization_mechanical_pass.f90"
#include "../homogenization_mechanical_isostrain.f90"
#include "../homogenization_mechanical_RGC.f90"
#include "../homogenization_thermal.f90"
#include "../homogenization_thermal_pass.f90"
#include "../homogenization_thermal_isotemperature.f90"
#include "../homogenization_damage.f90"
#include "../homogenization_damage_pass.f90"
#include "materialpoint_Marc.f90"

!--------------------------------------------------------------------------------------------------
!> @brief This is the MSC.Marc user subroutine for defining material behavior
!> @details (1) F,R,U are only available for continuum and membrane elements (not for
!> @details     shells and beams).
!> @details
!> @details (2) Use the -> 'Plasticity,3' card(=update+finite+large disp+constant d)
!> @details     in the parameter section of input deck (updated Lagrangian formulation).
!--------------------------------------------------------------------------------------------------
subroutine hypela2(d,g,e,de,s,t,dt,ngens,m,nn,kcus,matus,ndi,nshear,disp, &
                   dispt,coord,ffn,frotn,strechn,eigvn,ffn1,frotn1, &
                   strechn1,eigvn1,ncrd,itel,ndeg,ndm,nnode, &
                   jtype,lclass,ifr,ifu)
  use prec
  use DAMASK_interface
  use config
  use types
  use discretization_Marc
  use homogenization
  use materialpoint_Marc
  use OMP_LIB

  implicit none(type,external)
  integer(pI64),                         intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    ngens, &                                                                                        !< size of stress-strain law
    nn, &                                                                                           !< integration point number
    ndi, &                                                                                          !< number of direct components
    nshear, &                                                                                       !< number of shear components
    ncrd, &                                                                                         !< number of coordinates
    itel, &                                                                                         !< dimension of F and R, either 2 or 3
    ndeg, &                                                                                         !< number of degrees of freedom
    ndm, &                                                                                          !< not specified in MSC.Marc 2012 Manual D
    nnode, &                                                                                        !< number of nodes per element
    jtype, &                                                                                        !< element type
    ifr, &                                                                                          !< set to 1 if R has been calculated
    ifu                                                                                             !< set to 1 if stretch has been calculated
  integer(pI64), dimension(2),           intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    m, &                                                                                            !< (1) user element number, (2) internal element number
    matus, &                                                                                        !< (1) user material identification number, (2) internal material identification number
    kcus, &                                                                                         !< (1) layer number, (2) internal layer number
    lclass                                                                                          !< (1) element class, (2) 0: displacement, 1: low order Herrmann, 2: high order Herrmann
  real(pREAL),   dimension(*),           intent(in) :: &                                            ! has dimension(1) according to MSC.Marc 2012 Manual D, but according to example hypela2.f dimension(*)
    e, &                                                                                            !< total elastic strain
    de, &                                                                                           !< increment of strain
    dt                                                                                              !< increment of state variables
  real(pREAL),   dimension(itel),        intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    strechn, &                                                                                      !< square of principal stretch ratios, lambda(i) at t=n
    strechn1                                                                                        !< square of principal stretch ratios, lambda(i) at t=n+1
  real(pREAL),   dimension(3,3),         intent(in) :: &                                            ! has dimension(itel,*) according to MSC.Marc 2012 Manual D, but we alway assume dimension(3,3)
    ffn, &                                                                                          !< deformation gradient at t=n
    ffn1                                                                                            !< deformation gradient at t=n+1
  real(pREAL),   dimension(itel,*),      intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    frotn, &                                                                                        !< rotation tensor at t=n
    eigvn, &                                                                                        !< i principal direction components for j eigenvalues at t=n
    frotn1, &                                                                                       !< rotation tensor at t=n+1
    eigvn1                                                                                          !< i principal direction components for j eigenvalues at t=n+1
  real(pREAL),   dimension(ndeg,*),      intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    disp, &                                                                                         !< incremental displacements
    dispt                                                                                           !< displacements at t=n (at assembly, lovl=4) and displacements at t=n+1 (at stress recovery, lovl=6)
  real(pREAL),   dimension(ncrd,*),      intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    coord                                                                                           !< coordinates
  real(pREAL),   dimension(*),           intent(inout) :: &                                         ! according to MSC.Marc 2012 Manual D
    t                                                                                               !< state variables (comes in at t=n, must be updated to have state variables at t=n+1)
  real(pREAL),   dimension(ndi+nshear),  intent(out) :: &                                           ! has dimension(*) according to MSC.Marc 2012 Manual D, but we need to loop over it
    s, &                                                                                            !< stress - should be updated by user
    g                                                                                               !< change in stress due to temperature effects
  real(pREAL),   dimension(ngens,ngens), intent(out) :: &                                           ! according to MSC.Marc 2012 Manual D, but according to example hypela2.f dimension(ngens,*)
    d                                                                                               !< stress-strain law to be formed

!--------------------------------------------------------------------------------------------------
! Marc common blocks are in fixed format so they have to be reformated to free format (f90)
! Beware of changes in newer Marc versions

#include QUOTE(PASTE(include/concom,MARC4DAMASK))                                                   ! concom is needed for inc, lovl
#include QUOTE(PASTE(include/creeps,MARC4DAMASK))                                                   ! creeps is needed for timinc (time increment)

  logical :: cutBack
  real(pREAL), dimension(6) ::   stress
  real(pREAL), dimension(6,6) :: ddsdde
  integer :: computationMode, i, node, CPnodeID
  integer(pI32) :: defaultNumThreadsInt                                                             !< default value set by Marc

  integer, save :: &
    theInc       = -1, &                                                                            !< needs description
    lastLovl     =  0                                                                               !< lovl in previous call to marc hypela2
  real(pREAL), save :: &
    theTime      = 0.0_pREAL, &                                                                     !< needs description
    theDelta     = 0.0_pREAL
  logical, save :: &
    lastIncConverged  = .false., &                                                                  !< needs description
    outdatedByNewInc  = .false., &                                                                  !< needs description
    materialpoint_init_done   = .false.                                                             !< remember whether init has been done already


  defaultNumThreadsInt = omp_get_num_threads()                                                      ! remember number of threads set by Marc
  call omp_set_num_threads(1_pI32)                                                                  ! no openMP

  if (.not. materialpoint_init_done) then
    materialpoint_init_done = .true.
    call materialpoint_initAll()
  end if

  computationMode = 0                                                                               ! save initialization value, since it does not result in any calculation
  if (lovl == 4 ) then                                                                              ! jacobian requested by marc
    if (timinc < theDelta .and. theInc == inc .and. lastLovl /= lovl) &                             ! first after cutback
      computationMode = materialpoint_RESTOREJACOBIAN
  elseif (lovl == 6) then                                                                           ! stress requested by marc
    computationMode = materialpoint_CALCRESULTS
    if (cptim > theTime .or. inc /= theInc) then                                                    ! reached "convergence"
      cycleCounter = -1                                                                             ! first calc step increments this to cycle = 0
      if (inc == 0) then                                                                            ! >> start of analysis <<
        lastIncConverged = .false.
        outdatedByNewInc = .false.
        lastLovl = lovl                                                                             ! pretend that this is NOT the first after a lovl change
        print'(a,i6,1x,i2)', '<< HYPELA2 >> start of analysis..! ',m(1),nn
      else if (inc - theInc > 1) then                                                               ! >> restart of broken analysis <<
        lastIncConverged = .false.
        outdatedByNewInc = .false.
        print'(a,i6,1x,i2)', '<< HYPELA2 >> restart of analysis..! ',m(1),nn
      else                                                                                          ! >> just the next inc <<
        lastIncConverged = .true.
        outdatedByNewInc = .true.
        print'(a,i6,1x,i2)', '<< HYPELA2 >> new increment..! ',m(1),nn
      end if
    else if ( timinc < theDelta ) then                                                              ! >> cutBack <<
      lastIncConverged = .false.
      outdatedByNewInc = .false.
      cycleCounter = -1                                                                             ! first calc step increments this to cycle = 0
      print'(a,i6,1x,i2)', '<< HYPELA2 >> cutback detected..! ',m(1),nn
    end if                                                                                          ! convergence treatment end
    flush(6)

    if (lastLovl /= lovl) then
      cycleCounter  = cycleCounter + 1
      !mesh_cellnode = mesh_build_cellnodes()                                                       ! update cell node coordinates
      !call mesh_build_ipCoordinates()                                                              ! update ip coordinates
    end if
    if (outdatedByNewInc) then
      computationMode = ior(computationMode,materialpoint_AGERESULTS)
      outdatedByNewInc = .false.
    end if
    if (lastIncConverged) then
      computationMode = ior(computationMode,materialpoint_BACKUPJACOBIAN)
      lastIncConverged = .false.
    end if

    theTime  = cptim
    theDelta = timinc
    theInc   = inc

  end if
  lastLovl = lovl

  call materialpoint_general(computationMode,ffn,ffn1,t(1),timinc,int(m(1)),int(nn),stress,ddsdde)

  d = ddsdde(1:ngens,1:ngens)
  s = stress(1:ndi+nshear)
  g = 0.0_pREAL
  if (symmetricSolver) d = 0.5_pREAL*(d+transpose(d))

  call omp_set_num_threads(defaultNumThreadsInt)                                                    ! reset number of threads to stored default value

end subroutine hypela2


!--------------------------------------------------------------------------------------------------
!> @brief calculate internal heat generated due to inelastic energy dissipation
!--------------------------------------------------------------------------------------------------
subroutine flux(f,ts,n,time)
  use prec
  use homogenization
  use discretization_Marc

  implicit none(type,external)
  real(pREAL),   dimension(6),  intent(in) :: &
    ts
  integer(pI64), dimension(10), intent(in) :: &
    n
  real(pREAL),                  intent(in) :: &
    time
  real(pREAL),   dimension(2),  intent(out) :: &
    f


  f(1) = homogenization_f_T(discretization_Marc_FEM2DAMASK_cell(int(n(3)),int(n(1))))
  f(2) = 0.0_pREAL

 end subroutine flux


!--------------------------------------------------------------------------------------------------
!> @brief trigger writing of results
!> @details uedinc is called before each new increment, not at the end of a converged one.
!> Therefore, storing the last written inc with an 'save' variable is required to avoid writing the
! same increment multiple times.
!--------------------------------------------------------------------------------------------------
subroutine uedinc(inc,incsub)
  use prec
  use materialpoint_Marc
  use discretization_Marc

  implicit none(type,external)

  external :: nodvar
  integer(pI64), intent(in) :: inc, incsub

  integer :: n, nqncomp, nqdatatype
  integer, save :: inc_written
  real(pREAL), allocatable, dimension(:,:) :: d_n
#include QUOTE(PASTE(include/creeps,MARC4DAMASK))                                                   ! creeps is needed for timinc (time increment)


  if (inc > inc_written) then
    allocate(d_n(3,count(discretization_Marc_FEM2DAMASK_node /= -1)))
    do n = lbound(discretization_Marc_FEM2DAMASK_node,1), ubound(discretization_Marc_FEM2DAMASK_node,1)
      if (discretization_Marc_FEM2DAMASK_node(n) /= -1) then
        call nodvar(1,n,d_n(1:3,discretization_Marc_FEM2DAMASK_node(n)),nqncomp,nqdatatype)
        if (nqncomp == 2) d_n(3,discretization_Marc_FEM2DAMASK_node(n)) = 0.0_pREAL
      end if
    end do

    call discretization_Marc_UpdateNodeAndIpCoords(d_n)
    call materialpoint_result(int(inc),cptim)

    inc_written = int(inc)
  end if

end subroutine uedinc
