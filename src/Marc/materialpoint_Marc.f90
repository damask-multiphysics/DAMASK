!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief materialpoint engine
!--------------------------------------------------------------------------------------------------
module materialpoint_Marc
  use DAMASK_interface
  use prec
  use IO
  use YAML_types
  use YAML_parse
  use HDF5_utilities
  use results
  use config
  use math
  use rotations
  use polynomials
  use lattice
  use material
  use phase
  use homogenization

  use discretization
  use discretization_Marc

  implicit none
  private

  real(pReal), dimension (:,:,:),   allocatable, private :: &
    materialpoint_cs                                                                                !< Cauchy stress
  real(pReal), dimension (:,:,:,:), allocatable, private :: &
    materialpoint_dcsdE                                                                             !< Cauchy stress tangent
  real(pReal), dimension (:,:,:,:), allocatable, private :: &
    materialpoint_dcsdE_knownGood                                                                   !< known good tangent

  integer,                                       public :: &
    cycleCounter = 0                                                                                !< needs description

  integer, parameter,                            public :: &
    materialpoint_CALCRESULTS     = 2**0, &
    materialpoint_AGERESULTS      = 2**1, &
    materialpoint_BACKUPJACOBIAN  = 2**2, &
    materialpoint_RESTOREJACOBIAN = 2**3

  type, private :: tNumerics
    integer :: &
      iJacoStiffness                                                                                !< frequency of stiffness update
  end type tNumerics

  type(tNumerics), private :: num

  type, private :: tDebugOptions
    logical :: &
      basic, &
      extensive, &
      selective
    integer:: &
      element, &
      ip
  end type tDebugOptions

  type(tDebugOptions), private :: debugmaterialpoint

  public :: &
    materialpoint_general, &
    materialpoint_initAll, &
    materialpoint_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize all modules.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_initAll

  call DAMASK_interface_init
  call prec_init
  call IO_init
  call YAML_types_init
  call YAML_parse_init
  call HDF5_utilities_init
  call results_init(.false.)
  call config_init
  call math_init
  call rotations_init
  call polynomials_init
  call lattice_init
  call discretization_Marc_init
  call material_init(.false.)
  call phase_init
  call homogenization_init
  call materialpoint_init
  call config_deallocate

end subroutine materialpoint_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module materialpoint and initialize them
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_init

  class(tNode), pointer :: &
    debug_materialpoint

  print'(/,1x,a)', '<<<+-  materialpoint init  -+>>>'; flush(IO_STDOUT)

  allocate(materialpoint_cs(               6,discretization_nIPs,discretization_Nelems), source= 0.0_pReal)
  allocate(materialpoint_dcsdE(          6,6,discretization_nIPs,discretization_Nelems), source= 0.0_pReal)
  allocate(materialpoint_dcsdE_knownGood(6,6,discretization_nIPs,discretization_Nelems), source= 0.0_pReal)

!------------------------------------------------------------------------------
! read debug options

  debug_materialpoint => config_debug%get('materialpoint',defaultVal=emptyList)
  debugmaterialpoint%basic     = debug_materialpoint%contains('basic')
  debugmaterialpoint%extensive = debug_materialpoint%contains('extensive')
  debugmaterialpoint%selective = debug_materialpoint%contains('selective')
  debugmaterialpoint%element   = config_debug%get_asInt('element',defaultVal = 1)
  debugmaterialpoint%ip        = config_debug%get_asInt('integrationpoint',defaultVal = 1)

  if(debugmaterialpoint%basic) then
    print'(a32,1x,6(i8,1x))',   'materialpoint_cs:              ', shape(materialpoint_cs)
    print'(a32,1x,6(i8,1x))',   'materialpoint_dcsdE:           ', shape(materialpoint_dcsdE)
    print'(a32,1x,6(i8,1x),/)', 'materialpoint_dcsdE_knownGood: ', shape(materialpoint_dcsdE_knownGood)
    flush(IO_STDOUT)
  endif

end subroutine materialpoint_init


!--------------------------------------------------------------------------------------------------
!> @brief Update variables and call the material model.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_general(mode, ffn, ffn1, temperature_inp, dt, elFE, ip, cauchyStress, jacobian)

  integer, intent(in) ::                              elFE, &                                       !< FE element number
                                                      ip                                            !< integration point number
  real(pReal), intent(in) ::                          dt                                            !< time increment
  real(pReal), dimension (3,3), intent(in) ::         ffn, &                                        !< deformation gradient for t=t0
                                                      ffn1                                          !< deformation gradient for t=t1
  integer, intent(in) ::                              mode                                          !< computation mode  1: regular computation plus aging of results
  real(pReal), intent(in) ::                          temperature_inp                               !< temperature
  real(pReal), dimension(6), intent(out) ::           cauchyStress                                  !< stress as 6 vector
  real(pReal), dimension(6,6), intent(out) ::         jacobian                                      !< jacobian as 66 tensor (Consistent tangent dcs/dE)

  real(pReal)                                         J_inverse, &                                  ! inverse of Jacobian
                                                      rnd
  real(pReal), dimension (3,3) ::                     Kirchhoff                                     ! Piola-Kirchhoff stress
  real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                      H

  integer                                             elCP, &                                       ! crystal plasticity element number
                                                      i, j, k, l, m, n, ph, homog, mySource,ce

  real(pReal), parameter ::                          ODD_STRESS    = 1e15_pReal, &                  !< return value for stress if terminallyIll
                                                     ODD_JACOBIAN  = 1e50_pReal                     !< return value for jacobian if terminallyIll


  elCP = discretization_Marc_FEM2DAMASK_elem(elFE)
  ce   = discretization_Marc_FEM2DAMASK_cell(ip,elFE)

  if (debugmaterialpoint%basic .and. elCP == debugmaterialpoint%element .and. ip == debugmaterialpoint%ip) then
    print'(/,a)', '#############################################'
    print'(a1,a22,1x,i8,a13)',   '#','element',        elCP,         '#'
    print'(a1,a22,1x,i8,a13)',   '#','ip',             ip,           '#'
    print'(a1,a22,1x,i8,a13)',   '#','cycleCounter',   cycleCounter, '#'
    print'(a1,a22,1x,i8,a13)',   '#','computationMode',mode,         '#'
    if (terminallyIll) &
    print'(a,/)', '#           --- terminallyIll ---           #'
    print'(a,/)', '#############################################'; flush (6)
  endif

  if (iand(mode, materialpoint_BACKUPJACOBIAN) /= 0) &
    materialpoint_dcsde_knownGood = materialpoint_dcsde
  if (iand(mode, materialpoint_RESTOREJACOBIAN) /= 0) &
    materialpoint_dcsde = materialpoint_dcsde_knownGood

  if (iand(mode, materialpoint_AGERESULTS) /= 0) call materialpoint_forward

    homogenization_F0(1:3,1:3,ce) = ffn
    homogenization_F(1:3,1:3,ce) = ffn1

  if (iand(mode, materialpoint_CALCRESULTS) /= 0) then

    validCalculation: if (terminallyIll) then
      call random_number(rnd)
      if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
      materialpoint_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
      materialpoint_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_eye(6)

    else validCalculation
      if (debugmaterialpoint%extensive)  print'(a,i8,1x,i2)', '<< materialpoint >> calculation for elFE ip ',elFE,ip
      call homogenization_mechanical_response(dt,(elCP-1)*discretization_nIPs + ip,(elCP-1)*discretization_nIPs + ip)
      if (.not. terminallyIll) &
        call homogenization_mechanical_response2(dt,[ip,ip],[elCP,elCP])

      terminalIllness: if (terminallyIll) then

        call random_number(rnd)
        if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
        materialpoint_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
        materialpoint_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_eye(6)

      else terminalIllness

        ! translate from P to sigma
        Kirchhoff = matmul(homogenization_P(1:3,1:3,ce), transpose(homogenization_F(1:3,1:3,ce)))
        J_inverse  = 1.0_pReal / math_det33(homogenization_F(1:3,1:3,ce))
        materialpoint_cs(1:6,ip,elCP) = math_sym33to6(J_inverse * Kirchhoff,weighted=.false.)

        !  translate from dP/dF to dCS/dE
        H = 0.0_pReal
        do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
          H(i,j,k,l) = H(i,j,k,l) &
                     +  homogenization_F(j,m,ce) * homogenization_F(l,n,ce) &
                                                 * homogenization_dPdF(i,m,k,n,ce) &
                     -  math_delta(j,l) * homogenization_F(i,m,ce) * homogenization_P(k,m,ce) &
                     +  0.5_pReal * (  Kirchhoff(j,l)*math_delta(i,k) + Kirchhoff(i,k)*math_delta(j,l) &
                                     + Kirchhoff(j,k)*math_delta(i,l) + Kirchhoff(i,l)*math_delta(j,k))
        enddo; enddo; enddo; enddo; enddo; enddo

        forall(i=1:3, j=1:3,k=1:3,l=1:3) &
          H_sym(i,j,k,l) = 0.25_pReal * (H(i,j,k,l) + H(j,i,k,l) + H(i,j,l,k) + H(j,i,l,k))

        materialpoint_dcsde(1:6,1:6,ip,elCP) = math_sym3333to66(J_inverse * H_sym,weighted=.false.)

      endif terminalIllness
    endif validCalculation

    if (debugmaterialpoint%extensive &
        .and. ((debugmaterialpoint%element == elCP .and. debugmaterialpoint%ip == ip) .or. .not. debugmaterialpoint%selective)) then
        print'(a,i8,1x,i2,/,12x,6(f10.3,1x)/)', &
          '<< materialpoint >> stress/MPa at elFE ip ',   elFE, ip, materialpoint_cs(1:6,ip,elCP)*1.0e-6_pReal
        print'(a,i8,1x,i2,/,6(12x,6(f10.3,1x)/))', &
          '<< materialpoint >> Jacobian/GPa at elFE ip ', elFE, ip, transpose(materialpoint_dcsdE(1:6,1:6,ip,elCP))*1.0e-9_pReal
        flush(IO_STDOUT)
    endif

  endif

  if (all(abs(materialpoint_dcsdE(1:6,1:6,ip,elCP)) < 1e-10_pReal)) &
    call IO_warning(601,label1='element (CP)',ID1=elCP,label2='IP',ID2=ip)

  cauchyStress = materialpoint_cs   (1:6,    ip,elCP)
  jacobian     = materialpoint_dcsdE(1:6,1:6,ip,elCP)

end subroutine materialpoint_general


!--------------------------------------------------------------------------------------------------
!> @brief Forward data for new time increment.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_forward

  call homogenization_forward
  call phase_forward

end subroutine materialpoint_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_results(inc,time)

  integer,     intent(in) :: inc
  real(pReal), intent(in) :: time

  call results_openJobFile
  call results_addIncrement(inc,time)
  call phase_results
  call homogenization_results
  call discretization_results
  call results_finalizeIncrement
  call results_closeJobFile

end subroutine materialpoint_results

end module materialpoint_Marc
