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
    real(pReal), dimension(:,:,:), allocatable     :: &
      P_d, &
      P_t, &
      P_n
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

  integer :: Ninstance,p,i
  character(len=pStringLen) :: extmsg = ''
  real(pReal), dimension(:,:), allocatable :: d,n,t

  write(6,'(/,a)') ' <<<+-  kinematics_'//KINEMATICS_slipplane_opening_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_kinematics == KINEMATICS_slipplane_opening_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(kinematics_slipplane_opening_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    kinematics_slipplane_opening_instance(p) = count(phase_kinematics(:,1:p) == kinematics_slipplane_opening_ID)
    if (all(phase_kinematics(:,p) /= KINEMATICS_slipplane_opening_ID)) cycle
    associate(prm => param(kinematics_slipplane_opening_instance(p)), &
             config => config_phase(p))

    prm%sdot0    = config%getFloat('anisoductile_sdot0')
    prm%n        = config%getFloat('anisoductile_ratesensitivity')
    prm%Nslip    = config%getInts('nslip')
    prm%totalNslip = sum(prm%Nslip)

    d = lattice_slip_direction (prm%Nslip,config%getString('lattice_structure'),&
                                config%getFloat('c/a',defaultVal=0.0_pReal))
    t = lattice_slip_transverse(prm%Nslip,config%getString('lattice_structure'),&
                                config%getFloat('c/a',defaultVal=0.0_pReal))
    n = lattice_slip_normal    (prm%Nslip,config%getString('lattice_structure'),&
                                config%getFloat('c/a',defaultVal=0.0_pReal))
    allocate(prm%P_d(3,3,size(d,2)),prm%P_t(3,3,size(t,2)),prm%P_n(3,3,size(n,2)))

    do i=1, size(n,2)
      prm%P_d(1:3,1:3,i) = math_outer(d(1:3,i), n(1:3,i))
      prm%P_t(1:3,1:3,i) = math_outer(t(1:3,i), n(1:3,i))
      prm%P_n(1:3,1:3,i) = math_outer(n(1:3,i), n(1:3,i))
    enddo

    prm%critLoad = config%getFloats('anisoductile_criticalload',requiredSize=size(prm%Nslip))

    ! expand: family => system
    prm%critLoad = math_expand(prm%critLoad, prm%Nslip)

    ! sanity checks
    if (prm%n            <= 0.0_pReal)  extmsg = trim(extmsg)//' anisoDuctile_n'
    if (prm%sdot0        <= 0.0_pReal)  extmsg = trim(extmsg)//' anisoDuctile_sdot0'
    if (any(prm%critLoad <  0.0_pReal)) extmsg = trim(extmsg)//' anisoDuctile_critLoad'

!--------------------------------------------------------------------------------------------------
!  exit if any parameter is out of range
    if (extmsg /= '') call IO_error(211,ext_msg=trim(extmsg)//'('//SOURCE_DAMAGE_ANISODUCTILE_LABEL//')')

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

    traction_d = math_mul33xx33(S,prm%P_d(3,3,i))
    traction_t = math_mul33xx33(S,prm%P_t(3,3,i))
    traction_n = math_mul33xx33(S,prm%P_n(3,3,i))

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
      dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) &
                          + dudotd_dt*prm%P_d(k,l,i)*prm%P_d(m,n,i) &
                          + dudott_dt*prm%P_t(k,l,i)*prm%P_t(m,n,i) &
                          + dudotn_dt*prm%P_n(k,l,i)*prm%P_n(m,n,i)

    Ld = Ld &
       + udotd*prm%P_d(1:3,1:3,i) &
       + udott*prm%P_t(1:3,1:3,i) &
       + udotn*prm%P_n(1:3,1:3,i)
  enddo

  end associate

end subroutine kinematics_slipplane_opening_LiAndItsTangent

end module kinematics_slipplane_opening
