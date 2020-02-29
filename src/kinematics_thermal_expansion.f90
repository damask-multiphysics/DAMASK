!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from thermal expansion
!> @details to be done
!--------------------------------------------------------------------------------------------------
module kinematics_thermal_expansion
  use prec
  use IO
  use config
  use debug
  use math
  use lattice
  use material

  implicit none
  private

  integer, dimension(:), allocatable :: kinematics_thermal_expansion_instance

  type :: tParameters
    real(pReal) :: &
      T_ref
    real(pReal), allocatable, dimension(:,:,:) :: &
      expansion
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param

  public :: &
    kinematics_thermal_expansion_init, &
    kinematics_thermal_expansion_initialStrain, &
    kinematics_thermal_expansion_LiAndItsTangent

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_init

  integer :: Ninstance,p,i
  real(pReal),  dimension(:), allocatable :: &
    temp

  write(6,'(/,a)') ' <<<+-  kinematics_'//KINEMATICS_thermal_expansion_LABEL//' init  -+>>>'; flush(6)

  Ninstance = count(phase_kinematics == KINEMATICS_thermal_expansion_ID)
  if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
    write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

  allocate(kinematics_thermal_expansion_instance(size(config_phase)), source=0)
  allocate(param(Ninstance))

  do p = 1, size(config_phase)
    kinematics_thermal_expansion_instance(p) = count(phase_kinematics(:,1:p) == KINEMATICS_thermal_expansion_ID)
    if (all(phase_kinematics(:,p) /= KINEMATICS_thermal_expansion_ID)) cycle

    associate(prm => param(kinematics_thermal_expansion_instance(p)), &
              config => config_phase(p))

    prm%T_ref = config%getFloat('reference_temperature', defaultVal=0.0_pReal)
    ! read up to three parameters (constant, linear, quadratic with T)
    temp = config%getFloats('thermal_expansion11')
    !lattice_thermalExpansion33(1,1,1:size(temp),p) = temp
    temp = config%getFloats('thermal_expansion22', &
                            defaultVal=[(0.0_pReal, i=1,size(temp))],requiredSize=size(temp))
    !lattice_thermalExpansion33(2,2,1:size(temp),p) = temp
    temp = config%getFloats('thermal_expansion33', &
                            defaultVal=[(0.0_pReal, i=1,size(temp))],requiredSize=size(temp))
    end associate
  enddo

end subroutine kinematics_thermal_expansion_init


!--------------------------------------------------------------------------------------------------
!> @brief  report initial thermal strain based on current temperature deviation from reference
!--------------------------------------------------------------------------------------------------
pure function kinematics_thermal_expansion_initialStrain(homog,phase,offset)

  integer, intent(in) :: &
    phase, &
    homog, &
    offset
  real(pReal), dimension(3,3) :: &
    kinematics_thermal_expansion_initialStrain                                                      !< initial thermal strain (should be small strain, though)

  associate(prm => param(kinematics_thermal_expansion_instance(phase)))
  kinematics_thermal_expansion_initialStrain = &
    (temperature(homog)%p(offset) - prm%T_ref)**1 / 1. * &
    lattice_thermalExpansion33(1:3,1:3,1,phase) + &                                                 ! constant  coefficient
    (temperature(homog)%p(offset) - prm%T_ref)**2 / 2. * &
    lattice_thermalExpansion33(1:3,1:3,2,phase) + &                                                 ! linear    coefficient
    (temperature(homog)%p(offset) - prm%T_ref)**3 / 3. * &
    lattice_thermalExpansion33(1:3,1:3,3,phase)                                                     ! quadratic coefficient
  end associate

end function kinematics_thermal_expansion_initialStrain


!--------------------------------------------------------------------------------------------------
!> @brief  contains the constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
subroutine kinematics_thermal_expansion_LiAndItsTangent(Li, dLi_dTstar, ipc, ip, el)

  integer, intent(in) :: &
    ipc, &                                                                                          !< grain number
    ip, &                                                                                           !< integration point number
    el                                                                                              !< element number
  real(pReal),   intent(out), dimension(3,3) :: &
    Li                                                                                              !< thermal velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLi_dTstar                                                                                      !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)

  integer :: &
    phase, &
    homog
  real(pReal) :: &
    T, TDot

  phase = material_phaseAt(ipc,el)
  homog = material_homogenizationAt(el)
  T = temperature(homog)%p(thermalMapping(homog)%p(ip,el))
  TDot = temperatureRate(homog)%p(thermalMapping(homog)%p(ip,el))

  associate(prm => param(kinematics_thermal_expansion_instance(phase)))
  Li = TDot * ( &
                lattice_thermalExpansion33(1:3,1:3,1,phase)*(T - prm%T_ref)**0 &                    ! constant  coefficient
              + lattice_thermalExpansion33(1:3,1:3,2,phase)*(T - prm%T_ref)**1 &                    ! linear    coefficient
              + lattice_thermalExpansion33(1:3,1:3,3,phase)*(T - prm%T_ref)**2 &                    ! quadratic coefficient
              ) / &
       (1.0_pReal &
             + lattice_thermalExpansion33(1:3,1:3,1,phase)*(T - prm%T_ref)**1 / 1. &
             + lattice_thermalExpansion33(1:3,1:3,2,phase)*(T - prm%T_ref)**2 / 2. &
             + lattice_thermalExpansion33(1:3,1:3,3,phase)*(T - prm%T_ref)**3 / 3. &
       )
  end associate
  dLi_dTstar = 0.0_pReal

end subroutine kinematics_thermal_expansion_LiAndItsTangent

end module kinematics_thermal_expansion
