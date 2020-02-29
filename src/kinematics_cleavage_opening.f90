!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of cleavage planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_cleavage_opening
  use prec
  use IO
  use config
  use debug
  use math
  use lattice
  use material

  implicit none
  private

  integer, dimension(:), allocatable :: kinematics_cleavage_opening_instance

  type :: tParameters                                                                               !< container type for internal constitutive parameters
    integer :: &
      totalNcleavage
    integer, dimension(:),   allocatable :: &
      Ncleavage                                                                                     !< active number of cleavage systems per family
    real(pReal) :: &
      sdot0, &
      n
    real(pReal),   dimension(:),   allocatable :: &
      critDisp, &
      critLoad
    real(pReal), dimension(:,:,:,:), allocatable :: &
      cleavage_systems
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param                                             !< containers of constitutive parameters (len Ninstance)

  public :: &
    kinematics_cleavage_opening_init, &
    kinematics_cleavage_opening_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_cleavage_opening_init

  integer :: Ninstance,p,instance

  write(6,'(/,a)') ' <<<+-  kinematics_'//KINEMATICS_cleavage_opening_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_kinematics == KINEMATICS_cleavage_opening_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(kinematics_cleavage_opening_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    kinematics_cleavage_opening_instance(p) = count(phase_kinematics(:,1:p) == kinematics_cleavage_opening_ID)
    if (all(phase_kinematics(:,p) /= KINEMATICS_cleavage_opening_ID)) cycle

    associate(prm => param(kinematics_cleavage_opening_instance(p)), &
              config => config_phase(p))

    instance = kinematics_cleavage_opening_instance(p)

    prm%sdot0      = config%getFloat('anisobrittle_sdot0')
    prm%n   = config%getFloat('anisobrittle_ratesensitivity')

    prm%Ncleavage = config%getInts('ncleavage',defaultVal=emptyIntArray)

    prm%critDisp = config%getFloats('anisobrittle_criticaldisplacement',requiredSize=size(prm%Ncleavage))
    prm%critLoad = config%getFloats('anisobrittle_criticalload',requiredSize=size(prm%Ncleavage))


    prm%totalNcleavage = sum(prm%Ncleavage)
        prm%cleavage_systems  = lattice_SchmidMatrix_cleavage (prm%Ncleavage,config%getString('lattice_structure'),&
                                                     config%getFloat('c/a',defaultVal=0.0_pReal))

   ! if (kinematics_cleavage_opening_sdot_0(instance) <= 0.0_pReal) &
   !   call IO_error(211,el=instance,ext_msg='sdot_0 ('//KINEMATICS_cleavage_opening_LABEL//')')
   !   call IO_error(211,el=instance,ext_msg='critical_displacement ('//KINEMATICS_cleavage_opening_LABEL//')')
   ! if (any(kinematics_cleavage_opening_critLoad(1:size(tempInt),instance) < 0.0_pReal)) &
   !   call IO_error(211,el=instance,ext_msg='critical_load ('//KINEMATICS_cleavage_opening_LABEL//')')
   ! if (kinematics_cleavage_opening_N(instance) <= 0.0_pReal) &
   !   call IO_error(211,el=instance,ext_msg='rate_sensitivity ('//KINEMATICS_cleavage_opening_LABEL//')')

    end associate
  enddo

end subroutine kinematics_cleavage_opening_init

!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine kinematics_cleavage_opening_LiAndItsTangent(Ld, dLd_dTstar, S, ipc, ip, el)

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
  instance = kinematics_cleavage_opening_instance(phase)
  homog = material_homogenizationAt(el)
  damageOffset = damageMapping(homog)%p(ip,el)

  Ld = 0.0_pReal
  dLd_dTstar = 0.0_pReal
    do i = 1,param(instance)%totalNcleavage
      traction_d    = math_mul33xx33(S,param(instance)%cleavage_systems(1:3,1:3,1,i))
      traction_t    = math_mul33xx33(S,param(instance)%cleavage_systems(1:3,1:3,2,i))
      traction_n    = math_mul33xx33(S,param(instance)%cleavage_systems(1:3,1:3,3,i))
      traction_crit = param(instance)%critLoad(i)* &
                      damage(homog)%p(damageOffset)**2.0_pReal
      udotd = &
        sign(1.0_pReal,traction_d)* &
        param(instance)%sdot0* &
        (max(0.0_pReal, abs(traction_d) - traction_crit)/traction_crit)**param(instance)%n
      if (abs(udotd) > tol_math_check) then
        Ld = Ld + udotd*param(instance)%cleavage_systems(1:3,1:3,1,i)
        dudotd_dt = sign(1.0_pReal,traction_d)*udotd*param(instance)%n/ &
                    max(0.0_pReal, abs(traction_d) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
            dudotd_dt*param(instance)%cleavage_systems(k,l,1,i)* &
                      param(instance)%cleavage_systems(m,n,1,i)
      endif

      udott = &
        sign(1.0_pReal,traction_t)* &
        param(instance)%sdot0* &
        (max(0.0_pReal, abs(traction_t) - traction_crit)/traction_crit)**param(instance)%n
      if (abs(udott) > tol_math_check) then
        Ld = Ld + udott*param(instance)%cleavage_systems(1:3,1:3,2,i)
        dudott_dt = sign(1.0_pReal,traction_t)*udott*param(instance)%n/ &
                    max(0.0_pReal, abs(traction_t) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
            dudott_dt*param(instance)%cleavage_systems(k,l,2,i)* &
                      param(instance)%cleavage_systems(m,n,2,i)
      endif

      udotn = &
        sign(1.0_pReal,traction_n)* &
        param(instance)%sdot0* &
        (max(0.0_pReal, abs(traction_n) - traction_crit)/traction_crit)**param(instance)%n
      if (abs(udotn) > tol_math_check) then
        Ld = Ld + udotn*param(instance)%cleavage_systems(1:3,1:3,3,i)
        dudotn_dt = sign(1.0_pReal,traction_n)*udotn*param(instance)%n/ &
                    max(0.0_pReal, abs(traction_n) - traction_crit)
        forall (k=1:3,l=1:3,m=1:3,n=1:3) &
          dLd_dTstar(k,l,m,n) = dLd_dTstar(k,l,m,n) + &
            dudotn_dt*param(instance)%cleavage_systems(k,l,3,i)* &
                      param(instance)%cleavage_systems(m,n,3,i)
      endif
  enddo

end subroutine kinematics_cleavage_opening_LiAndItsTangent

end module kinematics_cleavage_opening
