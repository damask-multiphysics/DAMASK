submodule(phase:mechanics) eigendeformation

  interface
    module function kinematics_cleavage_opening_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_cleavage_opening_init

    module function kinematics_slipplane_opening_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_slipplane_opening_init

    module function kinematics_thermal_expansion_init(kinematics_length) result(myKinematics)
      integer, intent(in) :: kinematics_length
      logical, dimension(:,:), allocatable :: myKinematics
    end function kinematics_thermal_expansion_init
  end interface


contains


module subroutine eigendeformation_init(phases)

  class(tNode), pointer :: &
    phases

  integer :: &
    ph
  class(tNode), pointer :: &
    phase, &
   kinematics

  print'(/,a)', ' <<<+-  phase_mechanics_eigendeformation init  -+>>>'

!--------------------------------------------------------------------------------------------------
! initialize kinematic mechanisms
  allocate(phase_Nkinematics(phases%length),source = 0)
  do ph = 1,phases%length
    phase => phases%get(ph)
    kinematics => phase%get('kinematics',defaultVal=emptyList)
    phase_Nkinematics(ph) = kinematics%length
  enddo

  allocate(phase_kinematics(maxval(phase_Nkinematics),phases%length), source = KINEMATICS_undefined_ID)

  if(maxval(phase_Nkinematics) /= 0) then
    where(kinematics_cleavage_opening_init(maxval(phase_Nkinematics)))  phase_kinematics = KINEMATICS_cleavage_opening_ID
    where(kinematics_slipplane_opening_init(maxval(phase_Nkinematics))) phase_kinematics = KINEMATICS_slipplane_opening_ID
    where(kinematics_thermal_expansion_init(maxval(phase_Nkinematics))) phase_kinematics = KINEMATICS_thermal_expansion_ID
  endif

end subroutine eigendeformation_init


end submodule eigendeformation
