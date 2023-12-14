!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @brief materialpoint engine
!--------------------------------------------------------------------------------------------------
module materialpoint_Marc
  use DAMASK_interface
  use prec
  use IO
  use types
  use YAML
  use HDF5_utilities
  use result
  use config
  use math
  use rotations
  use polynomials
  use tables
  use crystal
  use material
  use phase
  use homogenization
  use constants

  use discretization
  use discretization_Marc

  implicit none(type,external)
  private

  real(pREAL), dimension (:,:,:),   allocatable :: &
    materialpoint_cs                                                                                !< Cauchy stress
  real(pREAL), dimension (:,:,:,:), allocatable :: &
    materialpoint_dcsdE, &                                                                          !< Cauchy stress tangent
    materialpoint_F                                                                                 !< deformation gradient
  real(pREAL), dimension (:,:,:,:), allocatable :: &
    materialpoint_dcsdE_knownGood                                                                   !< known good tangent

  integer,                                       public :: &
    cycleCounter = 0                                                                                !< needs description
  integer(kind(STATUS_OK)) :: &
    status

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

  public :: &
    materialpoint_general, &
    materialpoint_initAll, &
    materialpoint_result

contains


!--------------------------------------------------------------------------------------------------
!> @brief Initialize all modules.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_initAll()

  call DAMASK_interface_init()
  call prec_init()
  call IO_init()
  call types_init()
  call YAML_init()
  call HDF5_utilities_init()
  call result_init(.false.)
  call config_init()
  call math_init()
  call rotations_init()
  call polynomials_init()
  call tables_init()
  call crystal_init()
  call discretization_Marc_init()
  call material_init(.false.)
  call phase_init()
  call homogenization_init()
  call materialpoint_init()
  call config_material_deallocate()
  call config_numerics_deallocate()


end subroutine materialpoint_initAll


!--------------------------------------------------------------------------------------------------
!> @brief Allocate the arrays defined in module materialpoint and initialize them.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_init()

  print'(/,1x,a)', '<<<+-  materialpoint init  -+>>>'; flush(IO_STDOUT)

  allocate(materialpoint_F(              3,3,discretization_nIPs,discretization_Nelems), source= 0.0_pREAL)
  allocate(materialpoint_cs(               6,discretization_nIPs,discretization_Nelems), source= 0.0_pREAL)
  allocate(materialpoint_dcsdE(          6,6,discretization_nIPs,discretization_Nelems), source= 0.0_pREAL)
  allocate(materialpoint_dcsdE_knownGood(6,6,discretization_nIPs,discretization_Nelems), source= 0.0_pREAL)


end subroutine materialpoint_init


!--------------------------------------------------------------------------------------------------
!> @brief Update variables and call the material model.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_general(mode, ffn, ffn1, temperature_inp, dt, elFE, ip, cauchyStress, jacobian)

  integer, intent(in) ::                              elFE, &                                       !< FE element number
                                                      ip                                            !< integration point number
  real(pREAL), intent(in) ::                          dt                                            !< time increment
  real(pREAL), dimension (3,3), intent(in) ::         ffn, &                                        !< deformation gradient for t=t0
                                                      ffn1                                          !< deformation gradient for t=t1
  integer, intent(in) ::                              mode                                          !< computation mode  1: regular computation plus aging of results
  real(pREAL), intent(in) ::                          temperature_inp                               !< temperature
  real(pREAL), dimension(6), intent(out) ::           cauchyStress                                  !< stress as 6 vector
  real(pREAL), dimension(6,6), intent(out) ::         jacobian                                      !< jacobian as 66 tensor (Consistent tangent dcs/dE)

  real(pREAL)                                         J_inverse, &                                  ! inverse of Jacobian
                                                      rnd
  real(pREAL), dimension (3,3) ::                     Kirchhoff                                     ! Piola-Kirchhoff stress
  real(pREAL), dimension (3,3,3,3) ::                 H_sym, &
                                                      H

  integer                                             elCP, &                                       ! crystal plasticity element number
                                                      i, j, k, l, m, n, ph, homog, mySource,ce

  real(pREAL), parameter ::                          ODD_STRESS    = 1e15_pREAL, &                  !< return value for stress if broken
                                                     ODD_JACOBIAN  = 1e50_pREAL                     !< return value for jacobian if broken

  elCP = discretization_Marc_FEM2DAMASK_elem(elFE)
  ce   = discretization_Marc_FEM2DAMASK_cell(ip,elFE)


  if (iand(mode, materialpoint_BACKUPJACOBIAN) /= 0) &
    materialpoint_dcsde_knownGood = materialpoint_dcsde
  if (iand(mode, materialpoint_RESTOREJACOBIAN) /= 0) &
    materialpoint_dcsde = materialpoint_dcsde_knownGood

  if (iand(mode, materialpoint_AGERESULTS) /= 0) call materialpoint_forward()

    homogenization_F(1:3,1:3,ce) = ffn1
    materialpoint_F(1:3,1:3,ip,elCP) = ffn1

  if (iand(mode, materialpoint_CALCRESULTS) /= 0) then

    validCalculation: if (status /= STATUS_OK) then
      call random_number(rnd)
      if (rnd < 0.5_pREAL) rnd = rnd - 1.0_pREAL
      materialpoint_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
      materialpoint_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_eye(6)

    else validCalculation
      call homogenization_mechanical_response(status, &
                                              dt,(elCP-1)*discretization_nIPs + ip, (elCP-1)*discretization_nIPs + ip)

      terminalIllness: if (status /= STATUS_OK) then

        call random_number(rnd)
        if (rnd < 0.5_pREAL) rnd = rnd - 1.0_pREAL
        materialpoint_cs(1:6,ip,elCP)        = ODD_STRESS * rnd
        materialpoint_dcsde(1:6,1:6,ip,elCP) = ODD_JACOBIAN * math_eye(6)

      else terminalIllness

        ! translate from P to sigma
        Kirchhoff = matmul(homogenization_P(1:3,1:3,ce), transpose(materialpoint_F(1:3,1:3,ip,elCP)))
        J_inverse  = 1.0_pREAL / math_det33(materialpoint_F(1:3,1:3,ip,elCP))
        materialpoint_cs(1:6,ip,elCP) = math_sym33to6(J_inverse * Kirchhoff,weighted=.false.)

        !  translate from dP/dF to dCS/dE
        H = 0.0_pREAL
        do i=1,3; do j=1,3; do k=1,3; do l=1,3; do m=1,3; do n=1,3
          H(i,j,k,l) = H(i,j,k,l) &
                     +  materialpoint_F(j,m,ip,elCP) * materialpoint_F(l,n,ip,elCP) &
                                                     * homogenization_dPdF(i,m,k,n,ce) &
                     -  math_delta(j,l) * materialpoint_F(i,m,ip,elCP) * homogenization_P(k,m,ce) &
                     +  0.5_pREAL * (  Kirchhoff(j,l)*math_delta(i,k) + Kirchhoff(i,k)*math_delta(j,l) &
                                     + Kirchhoff(j,k)*math_delta(i,l) + Kirchhoff(i,l)*math_delta(j,k))
        end do; end do; end do; end do; end do; end do

        forall(i=1:3, j=1:3,k=1:3,l=1:3) &
          H_sym(i,j,k,l) = 0.25_pREAL * (H(i,j,k,l) + H(j,i,k,l) + H(i,j,l,k) + H(j,i,l,k))

        materialpoint_dcsde(1:6,1:6,ip,elCP) = math_sym3333to66(J_inverse * H_sym,weighted=.false.)

      end if terminalIllness
    end if validCalculation

  end if

  if (all(abs(materialpoint_dcsdE(1:6,1:6,ip,elCP)) < 1e-10_pREAL)) &
    call IO_warning(601,label1='element (CP)',ID1=elCP,label2='IP',ID2=ip)

  cauchyStress = materialpoint_cs   (1:6,    ip,elCP)
  jacobian     = materialpoint_dcsdE(1:6,1:6,ip,elCP)

end subroutine materialpoint_general


!--------------------------------------------------------------------------------------------------
!> @brief Forward data for new time increment.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_forward

  call homogenization_forward()
  call phase_forward()

end subroutine materialpoint_forward


!--------------------------------------------------------------------------------------------------
!> @brief Trigger writing of results.
!--------------------------------------------------------------------------------------------------
subroutine materialpoint_result(inc,time)

  integer,     intent(in) :: inc
  real(pREAL), intent(in) :: time

  call result_openJobFile()
  call result_addIncrement(inc,time)
  call phase_result()
  call homogenization_result()
  call discretization_result()
  call result_finalizeIncrement()
  call result_closeJobFile()

end subroutine materialpoint_result

end module materialpoint_Marc
