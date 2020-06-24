!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief CPFEM engine
!--------------------------------------------------------------------------------------------------
module CPFEM
  use prec
  use debug
  use FEsolving
  use math
  use rotations
  use YAML_types
  use YAML_parse
  use discretization_marc
  use material
  use config
  use crystallite
  use homogenization
  use IO
  use discretization
  use DAMASK_interface
  use numerics
  use HDF5_utilities
  use results
  use lattice
  use constitutive

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
      iJacoStiffness                                                                                  !< frequency of stiffness update
  end type tNumerics

  type(tNumerics), private ::  num
  
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
  call numerics_init
  call debug_init
  call config_init
  call math_init
  call rotations_init
  call YAML_types_init
  call YAML_init
  call HDF5_utilities_init
  call results_init(.false.)
  call discretization_marc_init
  call lattice_init
  call material_init(.false.)
  call constitutive_init
  call crystallite_init
  call homogenization_init
  call CPFEM_init

end subroutine CPFEM_initAll


!--------------------------------------------------------------------------------------------------
!> @brief allocate the arrays defined in module CPFEM and initialize them
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_init

  class(tNode), pointer :: &
    num_commercialFEM

  write(6,'(/,a)')   ' <<<+-  CPFEM init  -+>>>'
  flush(6)

  allocate(CPFEM_cs(               6,discretization_nIP,discretization_nElem), source= 0.0_pReal)
  allocate(CPFEM_dcsdE(          6,6,discretization_nIP,discretization_nElem), source= 0.0_pReal)
  allocate(CPFEM_dcsdE_knownGood(6,6,discretization_nIP,discretization_nElem), source= 0.0_pReal)

!------------------------------------------------------------------------------
! read numerical parameters and do sanity check 
  num_commercialFEM => numerics_root%get('commercialFEM',defaultVal=emptyDict)
  num%iJacoStiffness = num_commercialFEM%get_asInt('ijacostiffness',defaultVal=1)
  if (num%iJacoStiffness < 1)  call IO_error(301,ext_msg='iJacoStiffness')
!------------------------------------------------------------------------------

  if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0) then
    write(6,'(a32,1x,6(i8,1x))')   'CPFEM_cs:              ', shape(CPFEM_cs)
    write(6,'(a32,1x,6(i8,1x))')   'CPFEM_dcsdE:           ', shape(CPFEM_dcsdE)
    write(6,'(a32,1x,6(i8,1x),/)') 'CPFEM_dcsdE_knownGood: ', shape(CPFEM_dcsdE_knownGood)
    flush(6)
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
                                                      i, j, k, l, m, n, ph, homog, mySource
  logical                                             updateJaco                                    ! flag indicating if Jacobian has to be updated

  real(pReal), parameter ::                          ODD_STRESS    = 1e15_pReal, &                  !< return value for stress if terminallyIll
                                                     ODD_JACOBIAN  = 1e50_pReal                     !< return value for jacobian if terminallyIll


  elCP = mesh_FEM2DAMASK_elem(elFE)

  if (iand(debug_level(debug_CPFEM), debug_levelBasic) /= 0_pInt &
      .and. elCP == debug_e .and. ip == debug_i) then
    write(6,'(/,a)') '#############################################'
    write(6,'(a1,a22,1x,i8,a13)')   '#','element',        elCP,         '#'
    write(6,'(a1,a22,1x,i8,a13)')   '#','ip',             ip,           '#'
    write(6,'(a1,a22,1x,i8,a13)')   '#','cycleCounter',   cycleCounter, '#'
    write(6,'(a1,a22,1x,i8,a13)')   '#','computationMode',mode,         '#'
    if (terminallyIll) &
    write(6,'(a,/)') '#           --- terminallyIll ---           #'
    write(6,'(a,/)') '#############################################'; flush (6)
  endif

  if (iand(mode, CPFEM_BACKUPJACOBIAN) /= 0_pInt) &
    CPFEM_dcsde_knownGood = CPFEM_dcsde
  if (iand(mode, CPFEM_RESTOREJACOBIAN) /= 0_pInt) &
    CPFEM_dcsde = CPFEM_dcsde_knownGood

  if (iand(mode, CPFEM_AGERESULTS) /= 0_pInt) call CPFEM_forward

    chosenThermal1: select case (thermal_type(material_homogenizationAt(elCP)))
      case (THERMAL_conduction_ID) chosenThermal1
        temperature(material_homogenizationAt(elCP))%p(thermalMapping(material_homogenizationAt(elCP))%p(ip,elCP)) = &
          temperature_inp
      end select chosenThermal1
    materialpoint_F0(1:3,1:3,ip,elCP) = ffn
    materialpoint_F(1:3,1:3,ip,elCP) = ffn1

  if (iand(mode, CPFEM_CALCRESULTS) /= 0_pInt) then

    validCalculation: if (terminallyIll) then
      call random_number(rnd)
      if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
      CPFEM_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
      CPFEM_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_identity2nd(6)

    else validCalculation
      updateJaco = mod(cycleCounter,num%iJacoStiffness) == 0
      FEsolving_execElem = elCP
      FEsolving_execIP   = ip
      if (iand(debug_level(debug_CPFEM), debug_levelExtensive) /=  0_pInt) &
        write(6,'(a,i8,1x,i2)') '<< CPFEM >> calculation for elFE ip ',elFE,ip
      call materialpoint_stressAndItsTangent(updateJaco, dt)

      terminalIllness: if (terminallyIll) then

        call random_number(rnd)
        if (rnd < 0.5_pReal) rnd = rnd - 1.0_pReal
        CPFEM_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
        CPFEM_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_identity2nd(6)

      else terminalIllness

        ! translate from P to sigma
        Kirchhoff = matmul(materialpoint_P(1:3,1:3,ip,elCP), transpose(materialpoint_F(1:3,1:3,ip,elCP)))
        J_inverse  = 1.0_pReal / math_det33(materialpoint_F(1:3,1:3,ip,elCP))
        CPFEM_cs(1:6,ip,elCP) = math_sym33to6(J_inverse * Kirchhoff,weighted=.false.)

        !  translate from dP/dF to dCS/dE
        H = 0.0_pReal
        do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
          H(i,j,k,l) = H(i,j,k,l) &
                     +  materialpoint_F(j,m,ip,elCP) * materialpoint_F(l,n,ip,elCP) &
                                                     * materialpoint_dPdF(i,m,k,n,ip,elCP) &
                     -  math_delta(j,l) * materialpoint_F(i,m,ip,elCP) * materialpoint_P(k,m,ip,elCP) &
                     +  0.5_pReal * (  Kirchhoff(j,l)*math_delta(i,k) + Kirchhoff(i,k)*math_delta(j,l) &
                                     + Kirchhoff(j,k)*math_delta(i,l) + Kirchhoff(i,l)*math_delta(j,k))
        enddo; enddo; enddo; enddo; enddo; enddo

        forall(i=1:3, j=1:3,k=1:3,l=1:3) &
          H_sym(i,j,k,l) = 0.25_pReal * (H(i,j,k,l) + H(j,i,k,l) + H(i,j,l,k) + H(j,i,l,k))

        CPFEM_dcsde(1:6,1:6,ip,elCP) = math_sym3333to66(J_inverse * H_sym,weighted=.false.)

      endif terminalIllness
    endif validCalculation

    if ((iand(debug_level(debug_CPFEM), debug_levelExtensive) /= 0_pInt) &
         .and. ((debug_e == elCP .and. debug_i == ip) &
                .or. .not. iand(debug_level(debug_CPFEM), debug_levelSelective) /= 0_pInt)) then
        write(6,'(a,i8,1x,i2,/,12x,6(f10.3,1x)/)') &
          '<< CPFEM >> stress/MPa at elFE ip ',   elFE, ip, CPFEM_cs(1:6,ip,elCP)*1.0e-6_pReal
        write(6,'(a,i8,1x,i2,/,6(12x,6(f10.3,1x)/))') &
          '<< CPFEM >> Jacobian/GPa at elFE ip ', elFE, ip, transpose(CPFEM_dcsdE(1:6,1:6,ip,elCP))*1.0e-9_pReal
        flush(6)
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

  call crystallite_forward

end subroutine CPFEM_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine CPFEM_results(inc,time)

  integer(pInt), intent(in) :: inc
  real(pReal),   intent(in) :: time

  call results_openJobFile
  call results_addIncrement(inc,time)
  call constitutive_results
  call crystallite_results
  call homogenization_results
  call discretization_results
  call results_finalizeIncrement
  call results_closeJobFile

end subroutine CPFEM_results

end module CPFEM
