!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from thermal expansion
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:eigen) thermalexpansion

  integer, dimension(:), allocatable :: kinematics_thermal_expansion_instance

  type :: tParameters
    real(pReal) :: &
      T_ref
    real(pReal), dimension(3,3,3) :: &
      A = 0.0_pReal
  end type tParameters

  type(tParameters), dimension(:), allocatable :: param

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function thermalexpansion_init(kinematics_length) result(myKinematics)

  integer, intent(in)                  :: kinematics_length
  logical, dimension(:,:), allocatable :: myKinematics

  integer :: Ninstances,p,i,k
  real(pReal), dimension(:), allocatable :: temp
  class(tNode), pointer :: &
    phases, &
    phase, &
    mech, &
    kinematics, &
    kinematic_type

  print'(/,a)', ' <<<+-  phase:mechanical:eigen:thermalexpansion init  -+>>>'

  myKinematics = kinematics_active('thermalexpansion',kinematics_length)
  Ninstances = count(myKinematics)
  print'(a,i2)', ' # phases: ',Ninstances; flush(IO_STDOUT)
  if(Ninstances == 0) return

  phases => config_material%get('phase')
  allocate(param(Ninstances))
  allocate(kinematics_thermal_expansion_instance(phases%length), source=0)

  do p = 1, phases%length
    if(any(myKinematics(:,p))) kinematics_thermal_expansion_instance(p) = count(myKinematics(:,1:p))
    phase => phases%get(p)
    if(count(myKinematics(:,p)) == 0) cycle
    mech => phase%get('mechanical')
    kinematics => mech%get('eigen')
    do k = 1, kinematics%length
      if(myKinematics(k,p)) then
        associate(prm  => param(kinematics_thermal_expansion_instance(p)))
        kinematic_type => kinematics%get(k)

        prm%T_ref = kinematic_type%get_asFloat('T_ref', defaultVal=0.0_pReal)

        ! read up to three parameters (constant, linear, quadratic with T)
        temp = kinematic_type%get_as1dFloat('A_11')
        prm%A(1,1,1:size(temp)) = temp
        temp = kinematic_type%get_as1dFloat('A_33',defaultVal=[(0.0_pReal, i=1,size(temp))],requiredSize=size(temp))
        prm%A(3,3,1:size(temp)) = temp
        do i=1, size(prm%A,3)
          prm%A(1:3,1:3,i) = lattice_applyLatticeSymmetry33(prm%A(1:3,1:3,i),&
                                                   phase%get_asString('lattice'))
        enddo
        end associate
      endif
    enddo
  enddo


end function thermalexpansion_init


!--------------------------------------------------------------------------------------------------
!> @brief constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
module subroutine thermalexpansion_LiAndItsTangent(Li, dLi_dTstar, ph,me)

  integer, intent(in) :: ph, me
  real(pReal),   intent(out), dimension(3,3) :: &
    Li                                                                                              !< thermal velocity gradient
  real(pReal),   intent(out), dimension(3,3,3,3) :: &
    dLi_dTstar                                                                                      !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)

  real(pReal) :: T, dot_T

  
  T     = thermal_T(ph,me)
  dot_T = thermal_dot_T(ph,me)

  associate(prm => param(kinematics_thermal_expansion_instance(ph)))
    Li = dot_T * ( &
                  prm%A(1:3,1:3,1)*(T - prm%T_ref)**0 &                                             ! constant  coefficient
                + prm%A(1:3,1:3,2)*(T - prm%T_ref)**1 &                                             ! linear    coefficient
                + prm%A(1:3,1:3,3)*(T - prm%T_ref)**2 &                                             ! quadratic coefficient
                ) / &
         (1.0_pReal &
               + prm%A(1:3,1:3,1)*(T - prm%T_ref)**1 / 1. &
               + prm%A(1:3,1:3,2)*(T - prm%T_ref)**2 / 2. &
               + prm%A(1:3,1:3,3)*(T - prm%T_ref)**3 / 3. &
         )
  end associate
  dLi_dTstar = 0.0_pReal

end subroutine thermalexpansion_LiAndItsTangent

end submodule thermalexpansion
