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

#include "prec.f90"

module DAMASK_interface
  use prec
#if __INTEL_COMPILER >= 1800
  use, intrinsic :: ISO_fortran_env, only: &
    compiler_version, &
    compiler_options
#endif
  use ifport, only: &
    CHDIR

  implicit none
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

  print'(/,a)', ' <<<+-  DAMASK_marc init -+>>>'

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
  endif
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

  character(len=pStringLen) :: line
  integer :: myStat,fileUnit,s,e

  open(newunit=fileUnit, file=getSolverJobName()//INPUTFILEEXTENSION, &
       status='old', position='rewind', action='read',iostat=myStat)
  do
    read (fileUnit,'(A)',END=100) line
    if(index(trim(lc(line)),'solver') == 1) then
      read (fileUnit,'(A)',END=100) line                                                            ! next line
        s =     verify(line,      ' ')                                                              ! start of first chunk
        s = s + verify(line(s+1:),' ')                                                              ! start of second chunk
        e = s + scan  (line(s+1:),' ')                                                              ! end of second chunk
      solverIsSymmetric = line(s:e) /= '1'
    endif
  enddo
100 close(fileUnit)
  contains

  !--------------------------------------------------------------------------------------------------
  !> @brief changes characters in string to lower case
  !> @details copied from IO_lc
  !--------------------------------------------------------------------------------------------------
  function lc(string)

    character(len=*), intent(in) :: string                                                            !< string to convert
    character(len=len(string))   :: lc

    character(26), parameter :: LOWER = 'abcdefghijklmnopqrstuvwxyz'
    character(26), parameter :: UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer                  :: i,n

    do i=1,len(string)
      lc(i:i) = string(i:i)
      n = index(UPPER,lc(i:i))
      if (n/=0) lc(i:i) = LOWER(n:n)
    enddo
  end function lc

end function solverIsSymmetric

end module DAMASK_interface


#include "commercialFEM_fileList.f90"

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
  use YAML_types
  use discretization_marc
  use homogenization
  use CPFEM

  implicit none
  include "omp_lib.h"                                                                               ! the openMP function library
  integer,                               intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
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
  integer, dimension(2),                 intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    m, &                                                                                            !< (1) user element number, (2) internal element number
    matus, &                                                                                        !< (1) user material identification number, (2) internal material identification number
    kcus, &                                                                                         !< (1) layer number, (2) internal layer number
    lclass                                                                                          !< (1) element class, (2) 0: displacement, 1: low order Herrmann, 2: high order Herrmann
  real(pReal),   dimension(*),           intent(in) :: &                                            ! has dimension(1) according to MSC.Marc 2012 Manual D, but according to example hypela2.f dimension(*)
    e, &                                                                                            !< total elastic strain
    de, &                                                                                           !< increment of strain
    dt                                                                                              !< increment of state variables
  real(pReal),   dimension(itel),        intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    strechn, &                                                                                      !< square of principal stretch ratios, lambda(i) at t=n
    strechn1                                                                                        !< square of principal stretch ratios, lambda(i) at t=n+1
  real(pReal),   dimension(3,3),         intent(in) :: &                                            ! has dimension(itel,*) according to MSC.Marc 2012 Manual D, but we alway assume dimension(3,3)
    ffn, &                                                                                          !< deformation gradient at t=n
    ffn1                                                                                            !< deformation gradient at t=n+1
  real(pReal),   dimension(itel,*),      intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    frotn, &                                                                                        !< rotation tensor at t=n
    eigvn, &                                                                                        !< i principal direction components for j eigenvalues at t=n
    frotn1, &                                                                                       !< rotation tensor at t=n+1
    eigvn1                                                                                          !< i principal direction components for j eigenvalues at t=n+1
  real(pReal),   dimension(ndeg,*),      intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    disp, &                                                                                         !< incremental displacements
    dispt                                                                                           !< displacements at t=n (at assembly, lovl=4) and displacements at t=n+1 (at stress recovery, lovl=6)
  real(pReal),   dimension(ncrd,*),      intent(in) :: &                                            ! according to MSC.Marc 2012 Manual D
    coord                                                                                           !< coordinates
  real(pReal),   dimension(*),           intent(inout) :: &                                         ! according to MSC.Marc 2012 Manual D
    t                                                                                               !< state variables (comes in at t=n, must be updated to have state variables at t=n+1)
  real(pReal),   dimension(ndi+nshear),  intent(out) :: &                                           ! has dimension(*) according to MSC.Marc 2012 Manual D, but we need to loop over it
    s, &                                                                                            !< stress - should be updated by user
    g                                                                                               !< change in stress due to temperature effects
  real(pReal),   dimension(ngens,ngens), intent(out) :: &                                           ! according to MSC.Marc 2012 Manual D, but according to example hypela2.f dimension(ngens,*)
    d                                                                                               !< stress-strain law to be formed

!--------------------------------------------------------------------------------------------------
! Marc common blocks are in fixed format so they have to be reformated to free format (f90)
! Beware of changes in newer Marc versions

#include QUOTE(PASTE(./Marc/include/concom,Marc4DAMASK))                                            ! concom is needed for inc, lovl
#include QUOTE(PASTE(./Marc/include/creeps,Marc4DAMASK))                                            ! creeps is needed for timinc (time increment)

  logical :: cutBack
  real(pReal), dimension(6) ::   stress
  real(pReal), dimension(6,6) :: ddsdde
  integer :: computationMode, i, node, CPnodeID
  integer(pI32) :: defaultNumThreadsInt                                                             !< default value set by Marc

  integer(pInt), save :: &
    theInc       = -1_pInt, &                                                                       !< needs description
    lastLovl     =  0_pInt                                                                          !< lovl in previous call to marc hypela2
  real(pReal), save :: &
    theTime      = 0.0_pReal, &                                                                     !< needs description
    theDelta     = 0.0_pReal
  logical, save :: &
    lastIncConverged  = .false., &                                                                  !< needs description
    outdatedByNewInc  = .false., &                                                                  !< needs description
    CPFEM_init_done   = .false., &                                                                  !< remember whether init has been done already
    debug_basic       = .true.
  class(tNode), pointer :: &
    debug_Marc                                                                                      ! pointer to Marc debug options

  if(debug_basic) then
    print'(a,/,i8,i8,i2)', ' MSC.Marc information on shape of element(2), IP:', m, nn
    print'(a,2(i1))',      ' Jacobian:                      ', ngens,ngens
    print'(a,i1)',         ' Direct stress:                 ', ndi
    print'(a,i1)',         ' Shear stress:                  ', nshear
    print'(a,i2)',         ' DoF:                           ', ndeg
    print'(a,i2)',         ' Coordinates:                   ', ncrd
    print'(a,i12)',        ' Nodes:                         ', nnode
    print'(a,i1)',         ' Deformation gradient:          ', itel
    write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' Deformation gradient at t=n:', &
                                  transpose(ffn)
    write(6,'(/,a,/,3(3(f12.7,1x)/))',advance='no') ' Deformation gradient at t=n+1:', &
                                  transpose(ffn1)
  endif

  defaultNumThreadsInt = omp_get_num_threads()                                                      ! remember number of threads set by Marc
  call omp_set_num_threads(1_pI32)                                                                  ! no openMP

  if (.not. CPFEM_init_done) then
    CPFEM_init_done = .true.
    call CPFEM_initAll
    debug_Marc => config_debug%get('Marc',defaultVal=emptyList)
    debug_basic = debug_Marc%contains('basic')
  endif

  computationMode = 0                                                                               ! save initialization value, since it does not result in any calculation
  if (lovl == 4 ) then                                                                              ! jacobian requested by marc
    if (timinc < theDelta .and. theInc == inc .and. lastLovl /= lovl) &                             ! first after cutback
      computationMode = CPFEM_RESTOREJACOBIAN
  elseif (lovl == 6) then                                                                           ! stress requested by marc
    computationMode = CPFEM_CALCRESULTS
    if (cptim > theTime .or. inc /= theInc) then                                                    ! reached "convergence"
      terminallyIll = .false.
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
      endif
    else if ( timinc < theDelta ) then                                                              ! >> cutBack <<
      lastIncConverged = .false.
      outdatedByNewInc = .false.
      terminallyIll = .false.
      cycleCounter = -1                                                                             ! first calc step increments this to cycle = 0
      print'(a,i6,1x,i2)', '<< HYPELA2 >> cutback detected..! ',m(1),nn
    endif                                                                                           ! convergence treatment end
    flush(6)

    if (lastLovl /= lovl) then
      cycleCounter  = cycleCounter + 1
      !mesh_cellnode = mesh_build_cellnodes()                                                       ! update cell node coordinates
      !call mesh_build_ipCoordinates()                                                              ! update ip coordinates
    endif
    if (outdatedByNewInc) then
      computationMode = ior(computationMode,CPFEM_AGERESULTS)
      outdatedByNewInc = .false.
    endif
    if (lastIncConverged) then
      computationMode = ior(computationMode,CPFEM_BACKUPJACOBIAN)
      lastIncConverged = .false.
    endif

    theTime  = cptim
    theDelta = timinc
    theInc   = inc

  endif
  lastLovl = lovl

  call CPFEM_general(computationMode,ffn,ffn1,t(1),timinc,m(1),nn,stress,ddsdde)

  d = ddsdde(1:ngens,1:ngens)
  s = stress(1:ndi+nshear)
  g = 0.0_pReal
  if(symmetricSolver) d = 0.5_pReal*(d+transpose(d))

  call omp_set_num_threads(defaultNumThreadsInt)                                                    ! reset number of threads to stored default value

end subroutine hypela2


!--------------------------------------------------------------------------------------------------
!> @brief calculate internal heat generated due to inelastic energy dissipation
!--------------------------------------------------------------------------------------------------
subroutine flux(f,ts,n,time)
  use prec
  use homogenization
  use discretization_marc

  implicit none
  real(pReal), dimension(6),           intent(in) :: &
    ts
  integer,     dimension(10),          intent(in) :: &
    n
  real(pReal),                         intent(in) :: &
    time
  real(pReal), dimension(2),           intent(out) :: &
    f

  f(1) = homogenization_f_T(discretization_Marc_FEM2DAMASK_cell(n(3),n(1)))
  f(2) = 0.0_pReal

 end subroutine flux


!--------------------------------------------------------------------------------------------------
!> @brief trigger writing of results
!> @details uedinc is called before each new increment, not at the end of a converged one.
!> Therefore, storing the last written inc with an 'save' variable is required to avoid writing the
! same increment multiple times.
!--------------------------------------------------------------------------------------------------
subroutine uedinc(inc,incsub)
  use prec
  use CPFEM
  use discretization_marc

  implicit none
  integer, intent(in) :: inc, incsub
  integer :: n, nqncomp, nqdatatype
  integer, save :: inc_written
  real(pReal), allocatable, dimension(:,:) :: d_n
#include QUOTE(PASTE(./Marc/include/creeps,Marc4DAMASK))                                            ! creeps is needed for timinc (time increment)


  if (inc > inc_written) then
    allocate(d_n(3,count(discretization_Marc_FEM2DAMASK_node /= -1)))
    do n = lbound(discretization_Marc_FEM2DAMASK_node,1), ubound(discretization_Marc_FEM2DAMASK_node,1)
      if (discretization_Marc_FEM2DAMASK_node(n) /= -1) then
        call nodvar(1,n,d_n(1:3,discretization_Marc_FEM2DAMASK_node(n)),nqncomp,nqdatatype)
        if(nqncomp == 2) d_n(3,discretization_Marc_FEM2DAMASK_node(n)) = 0.0_pReal
      endif
    enddo

    call discretization_marc_UpdateNodeAndIpCoords(d_n)
    call CPFEM_results(inc,cptim)

    inc_written = inc
  endif

end subroutine uedinc
