submodule(constitutive) constitutive_damage

  implicit none

  interface

  module subroutine source_damage_anisoBrittle_init
  end subroutine source_damage_anisoBrittle_init 

  module subroutine source_damage_anisoDuctile_init
  end subroutine source_damage_anisoDuctile_init

  module subroutine source_damage_isoBrittle_init
  end subroutine source_damage_isoBrittle_init

  module subroutine source_damage_isoDuctile_init
  end subroutine source_damage_isoDuctile_init   

  module subroutine kinematics_cleavage_opening_init
  end subroutine kinematics_cleavage_opening_init

  module subroutine kinematics_slipplane_opening_init
  end subroutine kinematics_slipplane_opening_init


  module subroutine source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

    integer, intent(in) :: &
      phase, &
      constituent
    real(pReal),  intent(in) :: &
      phi
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi

  end subroutine source_damage_anisoBrittle_getRateAndItsTangent
 
  module subroutine source_damage_anisoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

    integer, intent(in) :: &
      phase, &
      constituent
    real(pReal),  intent(in) :: &
      phi
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi

  end subroutine source_damage_anisoDuctile_getRateAndItsTangent

  module subroutine source_damage_isoBrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

    integer, intent(in) :: &
      phase, &
      constituent
    real(pReal),  intent(in) :: &
      phi
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi

  end subroutine source_damage_isoBrittle_getRateAndItsTangent

  module subroutine source_damage_isoDuctile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

    integer, intent(in) :: &
      phase, &
      constituent
    real(pReal),  intent(in) :: &
      phi
    real(pReal),  intent(out) :: &
      localphiDot, &
      dLocalphiDot_dPhi

  end subroutine source_damage_isoDuctile_getRateAndItsTangent

  module subroutine source_thermal_dissipation_getRateAndItsTangent(TDot, dTDot_dT, Tstar, Lp, phase)

    integer, intent(in) :: &
      phase
    real(pReal),  intent(in), dimension(3,3) :: &
      Tstar
    real(pReal),  intent(in), dimension(3,3) :: &
      Lp

    real(pReal),  intent(out) :: &
      TDot, &
      dTDot_dT
  
  end subroutine source_thermal_dissipation_getRateAndItsTangent
 end interface

contains

module subroutine damage_init

! initialize source mechanisms
  if (any(phase_source == SOURCE_damage_isoBrittle_ID))       call source_damage_isoBrittle_init
  if (any(phase_source == SOURCE_damage_isoDuctile_ID))       call source_damage_isoDuctile_init
  if (any(phase_source == SOURCE_damage_anisoBrittle_ID))     call source_damage_anisoBrittle_init
  if (any(phase_source == SOURCE_damage_anisoDuctile_ID))     call source_damage_anisoDuctile_init

!--------------------------------------------------------------------------------------------------
! initialize kinematic mechanisms
  if (any(phase_kinematics == KINEMATICS_cleavage_opening_ID))  call kinematics_cleavage_opening_init
  if (any(phase_kinematics == KINEMATICS_slipplane_opening_ID)) call kinematics_slipplane_opening_init

end subroutine damage_init


module procedure constitutive_damage_getRateAndItsTangents

  real(pReal) :: &
    localphiDot, &
    dLocalphiDot_dPhi
  integer :: &
    phase, &
    grain, &
    source, &
    constituent

   phiDot = 0.0_pReal
   dPhiDot_dPhi = 0.0_pReal
 
   do grain = 1, homogenization_Ngrains(material_homogenizationAt(el))
     phase = material_phaseAt(grain,el)
     constituent = material_phasememberAt(grain,ip,el)
     do source = 1, phase_Nsources(phase)
       select case(phase_source(source,phase))
         case (SOURCE_damage_isoBrittle_ID)
           call source_damage_isobrittle_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (SOURCE_damage_isoDuctile_ID)
           call source_damage_isoductile_getRateAndItsTangent  (localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (SOURCE_damage_anisoBrittle_ID)
           call source_damage_anisobrittle_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case (SOURCE_damage_anisoDuctile_ID)
           call source_damage_anisoductile_getRateAndItsTangent(localphiDot, dLocalphiDot_dPhi, phi, phase, constituent)

         case default
         localphiDot = 0.0_pReal
         dLocalphiDot_dPhi = 0.0_pReal

      end select
      phiDot = phiDot + localphiDot
      dPhiDot_dPhi = dPhiDot_dPhi + dLocalphiDot_dPhi
    enddo
  enddo

end procedure constitutive_damage_getRateAndItsTangents

end submodule
