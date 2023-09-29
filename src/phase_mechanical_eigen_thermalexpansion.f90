!--------------------------------------------------------------------------------------------------
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from thermal expansion
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:eigen) thermalexpansion

  integer, dimension(:), allocatable :: kinematics_thermal_expansion_instance

  type :: tParameters
    type(tPolynomial) :: &
      Alpha_11, &
      Alpha_33
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

  integer :: p, k
  type(tList), pointer :: &
    kinematics
  type(tDict), pointer :: &
    phases, &
    phase, &
    mech


  myKinematics = kinematics_active('thermalexpansion',kinematics_length)
  if (count(myKinematics) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:eigen:thermalexpansion init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(myKinematics); flush(IO_STDOUT)

  phases => config_material%get_dict('phase')
  allocate(param(count(myKinematics)))
  allocate(kinematics_thermal_expansion_instance(phases%length), source=0)

  do p = 1, phases%length
    if (any(myKinematics(:,p))) kinematics_thermal_expansion_instance(p) = count(myKinematics(:,1:p))
    phase => phases%get_dict(p)
    if (count(myKinematics(:,p)) == 0) cycle
    mech => phase%get_dict('mechanical')
    kinematics => mech%get_list('eigen')
    do k = 1, kinematics%length
      if (myKinematics(k,p)) then
        associate(prm  => param(kinematics_thermal_expansion_instance(p)))

          prm%Alpha_11 = polynomial(kinematics%get_dict(k),'Alpha_11','T')
          if (any(phase_lattice(p) == ['hP','tI'])) &
            prm%Alpha_33 = polynomial(kinematics%get_dict(k),'Alpha_33','T')
        end associate
      end if
    end do
  end do

end function thermalexpansion_init


!--------------------------------------------------------------------------------------------------
!> @brief constitutive equation for calculating the velocity gradient
!--------------------------------------------------------------------------------------------------
module subroutine thermalexpansion_LiAndItsTangent(Li, dLi_dTstar, ph,me)

  integer, intent(in) :: ph, me
  real(pREAL),   intent(out), dimension(3,3) :: &
    Li                                                                                              !< thermal velocity gradient
  real(pREAL),   intent(out), dimension(3,3,3,3) :: &
    dLi_dTstar                                                                                      !< derivative of Li with respect to Tstar (4th-order tensor defined to be zero)

  real(pREAL) :: T, dot_T
  real(pREAL), dimension(3,3) :: Alpha


  T     = thermal_T(ph,me)
  dot_T = thermal_dot_T(ph,me)

  associate(prm => param(kinematics_thermal_expansion_instance(ph)))

    Alpha = 0.0_pREAL
    Alpha(1,1) = prm%Alpha_11%at(T)
    if (any(phase_lattice(ph) == ['hP','tI'])) Alpha(3,3) = prm%Alpha_33%at(T)
    Alpha = crystal_symmetrize_33(Alpha,phase_lattice(ph))
    Li = dot_T * Alpha

  end associate
  dLi_dTstar = 0.0_pREAL

end subroutine thermalexpansion_LiAndItsTangent

end submodule thermalexpansion
