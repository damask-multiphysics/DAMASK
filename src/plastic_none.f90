!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Dummy plasticity for purely elastic material
!--------------------------------------------------------------------------------------------------
submodule(constitutive) plastic_none

contains

!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module subroutine plastic_none_init

 integer :: &
   Ninstance, &
   p, &
   NipcMyPhase

 write(6,'(/,a)')   ' <<<+-  plastic_'//PLASTICITY_NONE_label//' init  -+>>>'

 Ninstance = count(phase_plasticity == PLASTICITY_NONE_ID)
 if (iand(debug_level(debug_constitutive),debug_levelBasic) /= 0) &
   write(6,'(a16,1x,i5,/)') '# instances:',Ninstance

 do p = 1, size(phase_plasticity)
   if (phase_plasticity(p) /= PLASTICITY_NONE_ID) cycle

   NipcMyPhase = count(material_phaseAt == p) * discretization_nIP
   call material_allocatePlasticState(p,NipcMyPhase,0,0,0, &
                                      0,0,0)
 enddo

end subroutine plastic_none_init

end submodule plastic_none
