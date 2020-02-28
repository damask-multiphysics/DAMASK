!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of slip planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_slipplane_opening
  use prec
  use config
  use IO
  use debug
  use math
  use lattice
  use material

  implicit none
  private

  integer, dimension(:), allocatable :: kinematics_slipplane_opening_instance

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      totalNslip
    integer, dimension(:),   allocatable :: &
      Nslip                                                                                         !< active number of slip systems per family
    real(pReal) :: &
      sdot0, &
      n
    real(pReal), dimension(:),   allocatable :: &
      critLoad
    real(pReal), dimension(:,:), allocatable     :: &
     slip_direction, &
     slip_normal, &
     slip_transverse
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstance)

  public :: &
    kinematics_slipplane_opening_init, &
    kinematics_slipplane_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_init

  integer :: Ninstance,p

  write(6,'(/,a)') ' <<<+-  kinematics_'//KINEMATICS_slipplane_opening_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_kinematics == KINEMATICS_slipplane_opening_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(kinematics_slipplane_opening_instance(size(config_phase)), source=0)
  do p = 1, size(config_phase)
    kinematics_slipplane_opening_instance(p) = count(phase_kinematics(:,1:p) == kinematics_slipplane_opening_ID) ! ToDo: count correct?
  enddo

  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    if (all(phase_kinematics(:,p) /= KINEMATICS_slipplane_opening_ID)) cycle
    associate(prm => param(kinematics_slipplane_opening_instance(p)), &
             config => config_phase(p))

    prm%sdot0    = config%getFloat('anisoductile_sdot0')
    prm%n        = config%getFloat('anisoductile_ratesensitivity')
    prm%Nslip    = config%getInts('nslip')

    prm%slip_direction  = lattice_slip_direction (prm%Nslip,config%getString('lattice_structure'),&
                                                  config%getFloat('c/a',defaultVal=0.0_pReal))
    prm%slip_normal     = lattice_slip_normal    (prm%Nslip,config%getString('lattice_structure'),&
                                                  config%getFloat('c/a',defaultVal=0.0_pReal))
    prm%slip_transverse = lattice_slip_transverse(prm%Nslip,config%getString('lattice_structure'),&
                                                  config%getFloat('c/a',defaultVal=0.0_pReal))

    prm%critLoad = config%getFloats('anisoductile_criticalload',requiredSize=size(prm%Nslip))
    prm%critLoad = math_expand(prm%critLoad, prm%Nslip)

    !  if (kinematics_slipplane_opening_sdot_0(instance) <= 0.0_pReal) &
    !    call IO_error(211,el=instance,ext_msg='sdot_0 ('//KINEMATICS_slipplane_opening_LABEL//')')
    !  if (any(kinematics_slipplane_opening_critPlasticStrain(:,instance) < 0.0_pReal)) &
    !    call IO_error(211,el=instance,ext_msg='criticaPlasticStrain ('//KINEMATICS_slipplane_opening_LABEL//')')
    !  if (kinematics_slipplane_opening_N(instance) <= 0.0_pReal) &
    !    call IO_error(211,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_slipplane_opening_LABEL//')')

    end associate
  enddo

end subroutine kinematics_slipplane_opening_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine kinematics_slipplane_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< grain number
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(in),  dimension(3,3) :: &
    S
  real(pReal),   intent(out), dimension(3,3) :: &
    Ld                                                                                              !< damage velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLd_dTstar                                                                                      !< derivative of Ld with respect to Tstar (4th-order tensor)
  real(pReal),   dimension(3,3) :: &
    projection_d, projection_t, projection_n                                                        !< projection modes 3x3 tensor
  integer :: &
    instance, phase, &
    homog, damageOffset, &
    i, k, l, m, n
  real(pReal) :: &
    traction_d, traction_t, traction_n, traction_crit, &
    udotd, dudotd_dt, udott, dudott_dt, udotn, dudotn_dt

  phase = material_phaseAt(ipc,el)
  instance = kinematics_slipplane_opening_instance(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)

  associate(prm => param(instance))
  Ld = 0.0_pReal
  dLd_dTstar = 0.0_pReal
  do i = 1, prm%totalNslip

    projection_d = math_outer(prm%slip_direction(1:3,i), prm%slip_normal(1:3,i))
    projection_t = math_outer(prm%slip_transverse(1:3,i),prm%slip_normal(1:3,i))
    projection_n = math_outer(prm%slip_normal(1:3,i),    prm%slip_normal(1:3,i))

    traction_d    = math_mul33xx33(S,projection_d)
    traction_t    = math_mul33xx33(S,projection_t)
    traction_n    = math_mul33xx33(S,projection_n)

    traction_crit = prm%critLoad(i)* damage(homog)%p(damageOffset)                                  ! degrading critical load carrying capacity by damage

    udotd = sign(1.0_pReal,traction_d)* prm%sdot0* (  abs(traction_d)/traction_crit &
                                                    - abs(traction_d)/prm%critLoad(i))**prm%n
    udott = sign(1.0_pReal,traction_t)* prm%sdot0* (  abs(traction_t)/traction_crit &
                                                    - abs(traction_t)/prm%critLoad(i))**prm%n
    udotn = prm%sdot0* (  max(0.0_pReal,traction_n)/traction_crit &
                        - max(0.0_pReal,traction_n)/prm%critLoad(i))**prm%n

    if (dNeq0(traction_d)) then
      dudotd_dt = udotd*prm%n/traction_d
    else
      dudotd_dt = 0.0_pReal
    endif
    if (dNeq0(traction_t)) then
      dudott_dt = udott*prm%n/traction_t
    else
      dudott_dt = 0.0_pReal
    endif
    if (dNeq0(traction_n)) then
      dudotn_dt = udotn*prm%n/traction_n
    else
      dudotn_dt = 0.0_pReal
    endif

    forall (k=1:3,l=1:3,m=1:3,n=1:3) &
      dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + dudotd_dt*projection_d(k,l)*projection_d(m,n) &
                                                + dudott_dt*projection_t(k,l)*projection_t(m,n) &
                                                + dudotn_dt*projection_n(k,l)*projection_n(m,n)

    Ld = Ld + udotd*projection_d &
            + udott*projection_t &
            + udotn*projection_n
  enddo

  end associate

end subroutine kinematics_slipplane_opening_LiAndItsTangent

end module kinematics_slipplane_opening
