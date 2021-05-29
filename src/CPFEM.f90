!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief CPFEM engine
!--------------------------------------------------------------------------------------------------
module CPFEM
  use prec
  use math
  use rotations
  use YAML_types
  use YAML_parse
  use discretization_marc
  use material
  use config
  use homogenization
  use IO
  use discretization
  use DAMASK_interface
  use HDF5_utilities
  use results
  use lattice
  use phase

  implicit none
  private

  real(pReal), dimension (:,:,:),   allocatable, private :: &
    CPFEM_cs                                                                                        !< Cauchy stress
  real(pReal), dimension (:,:,:,:), allocatable, private :: &
    CPFEM_dcsdE                                                                                     !< Cauchy stress tangent
  real(pReal), dimension (:,:,:,:), allocatable, private :: &
    CPFEM_dcsdE_knownGood                                                                           !< known good tangent

  integer(pInt),                                 public :: &
    cycleCounter =  0_pInt                                                                          !< needs description

  integer(pInt), parameter,                      public :: &
    CPFEM_CALCRESULTS     = 2_pInt**0_pInt, &
    CPFEM_AGERESULTS      = 2_pInt**1_pInt, &
    CPFEM_BACKUPJACOBIAN  = 2_pInt**2_pInt, &
    CPFEM_RESTOREJACOBIAN = 2_pInt**3_pInt

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

  type(tDebugOptions), private :: debugCPFEM

  public :: &
    CPFEM_general, &
    CPFEM_initAll, &
    CPFEM_results

contains


!--------------------------------------------------------------------------------------------------
!> @brief call all module initializations
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_initAll

  call DAMASK_interface_init
  call prec_init
  call IO_init
  call YAML_types_init
  call YAML_parse_init
  call config_init
  call math_init
  call rotations_init
  call HDF5_utilities_init
  call results_init(.false.)
  call discretization_marc_init
  call lattice_init
  call material_init(.false.)
  call phase_init
  call homogenization_init
  call crystallite_init
  call CPFEM_init
  call config_deallocate

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init

  class(tNode), pointer :: &
    debug_CPFEM

  print'(/,a)', ' <<<+-  CPFEM init  -+>>>'; flush(IO_STDOUT)

  allocate(CPFEM_cs(               6,discretization_nIPs,discretization_Nelems), source= 0.0_pReal)
  allocate(CPFEM_dcsdE(          6,6,discretization_nIPs,discretization_Nelems), source= 0.0_pReal)
  allocate(CPFEM_dcsdE_knownGood(6,6,discretization_nIPs,discretization_Nelems), source= 0.0_pReal)

!------------------------------------------------------------------------------
! read debug options

  debug_CPFEM => config_debug%get('CPFEM',defaultVal=emptyList)
  debugCPFEM%basic     = debug_CPFEM%contains('basic')
  debugCPFEM%extensive = debug_CPFEM%contains('extensive')
  debugCPFEM%selective = debug_CPFEM%contains('selective')
  debugCPFEM%element   = config_debug%get_asInt('element',defaultVal = 1)
  debugCPFEM%ip        = config_debug%get_asInt('integrationpoint',defaultVal = 1)

  if(debugCPFEM%basic) then
    print'(a32,1x,6(i8,1x))',   'CPFEM_cs:              ', shape(CPFEM_cs)
    print'(a32,1x,6(i8,1x))',   'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
    print'(a32,1x,6(i8,1x),/)', 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
    flush(IO_STDOUT)
  endif

end subroutine CPFEM_init


!--------------------------------------------------------------------------------------------------
!> @brief perform initialization at first call, update variables and call the actual material model
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_general(mode, ffn, ffn1, temperature_inp, dt, elFE, ip, cauchyStress, jacobian)

  integer(pInt), intent(in) ::                        elFE, &                                       !< FE element number
                                                      ip                                            !< integration point number
  real(pReal), intent(in) ::                          dt                                            !< time increment
  real(pReal), dimension (3,3), intent(in) ::         ffn, &                                        !< deformation gradient for t=t0
                                                      ffn1                                          !< deformation gradient for t=t1
  integer(pInt), intent(in) ::                        mode                                          !< computation mode  1: regular computation plus aging of results
  real(pReal), intent(in) ::                          temperature_inp                               !< temperature
  real(pReal), dimension(6), intent(out) ::           cauchyStress                                  !< stress as 6 vector
  real(pReal), dimension(6,6), intent(out) ::         jacobian                                      !< jacobian as 66 tensor (Consistent tangent dcs/dE)

  real(pReal)                                         J_inverse, &                                  ! inverse of Jacobian
                                                      rnd
  real(pReal), dimension (3,3) ::                     Kirchhoff                                     ! Piola-Kirchhoff stress
  real(pReal), dimension (3,3,3,3) ::                 H_sym, &
                                                      H

  integer(pInt)                                       elCP, &                                       ! crystal plasticity element number
                                                      i, j, k, l, m, n, ph, homog, mySource,ce

  real(pReal), parameter ::                          ODD_STRESS    = 1e15_pReal, &                  !< return value for stress if terminallyIll
                                                     ODD_JACOBIAN  = 1e50_pReal                     !< return value for jacobian if terminallyIll


  elCP = discretization_Marc_FEM2DAMASK_elem(elFE)
  ce   = discretization_Marc_FEM2DAMASK_cell(ip,elFE)

  if (debugCPFEM%basic .and. elCP == debugCPFEM%element .and. ip == debugCPFEM%ip) then
    print'(/,a)', '#############################################'
    print'(a1,a22,1x,i8,a13)',   '#','element',        elCP,         '#'
    print'(a1,a22,1x,i8,a13)',   '#','ip',             ip,           '#'
    print'(a1,a22,1x,i8,a13)',   '#','cycleCounter',   cycleCounter, '#'
    print'(a1,a22,1x,i8,a13)',   '#','computationMode',mode,         '#'
    if (terminallyIll) &
    print'(a,/)', '#           --- terminallyIll ---           #'
    print'(a,/)', '#############################################'; flush (6)
  endif

  if (iand(mode, CPFEM_BACKUPJACOBIAN) /= 0_pInt) &
    CPFEM_dcsde_knownGood = CPFEM_dcsde
  if (iand(mode, CPFEM_RESTOREJACOBIAN) /= 0_pInt) &
    CPFEM_dcsde = CPFEM_dcsde_knownGood

  if (iand(mode, CPFEM_AGERESULTS) /= 0_pInt) call CPFEM_forward

    homogenization_F0(1:3,1:3,ce) = ffn
    homogenization_F(1:3,1:3,ce) = ffn1

  if (iand(mode, CPFEM_CALCRESULTS) /= 0_pInt) then

    validCalculation: if (terminallyIll) then
      call random_number(rnd)
      if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
      CPFEM_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
      CPFEM_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_eye(6)

    else validCalculation
      if (debugCPFEM%extensive) &
        print'(a,i8,1x,i2)', '<< CPFEM >> calculation for elFE ip ',elFE,ip
      call materialpoint_stressAndItsTangent(dt,[ip,ip],[elCP,elCP])

      terminalIllness: if (terminallyIll) then

        call random_number(rnd)
        if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
        CPFEM_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
        CPFEM_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_eye(6)

      else terminalIllness

        ! translate from P to sigma
        Kirchhoff = matmul(homogenization_P(1:3,1:3,ce), transpose(homogenization_F(1:3,1:3,ce)))
        J_inverse  = 1.0_pReal / math_det33(homogenization_F(1:3,1:3,ce))
        CPFEM_cs(1:6,ip,elCP) = math_sym33to6(J_inverse * Kirchhoff,weighted=.false.)

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

        CPFEM_dcsde(1:6,1:6,ip,elCP) = math_sym3333to66(J_inverse * H_sym,weighted=.false.)

      endif terminalIllness
    endif validCalculation

    if (debugCPFEM%extensive &
         .and. ((debugCPFEM%element == elCP .and. debugCPFEM%ip == ip) .or. .not. debugCPFEM%selective)) then
        print'(a,i8,1x,i2,/,12x,6(f10.3,1x)/)', &
          '<< CPFEM >> stress/MPa at elFE ip ',   elFE, ip, CPFEM_cs(1:6,ip,elCP)*1.0e-6_pReal
        print'(a,i8,1x,i2,/,6(12x,6(f10.3,1x)/))', &
          '<< CPFEM >> Jacobian/GPa at elFE ip ', elFE, ip, transpose(CPFEM_dcsdE(1:6,1:6,ip,elCP))*1.0e-9_pReal
        flush(IO_STDOUT)
    endif

  endif

  if (all(abs(CPFEM_dcsdE(1:6,1:6,ip,elCP)) < 1e-10_pReal)) call IO_warning(601,elCP,ip)

  cauchyStress = CPFEM_cs   (1:6,    ip,elCP)
  jacobian     = CPFEM_dcsdE(1:6,1:6,ip,elCP)

end subroutine CPFEM_general


!--------------------------------------------------------------------------------------------------
!> @brief Forward data for new time increment.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_forward

  call homogenization_forward
  call phase_forward

end subroutine CPFEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)

  integer(pInt), intent(in) :: inc
  real(pReal),   intent(in) :: time

  call results_openJobFile
  call results_addIncrement(inc,time)
  call phase_results
  call homogenization_results
  call discretization_results
  call results_finalizeIncrement
  call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM
