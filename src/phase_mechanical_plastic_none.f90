!--------------------------------------------------------------------------------------------------
!> @author Franz Roters, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Max-Planck-Institut für Eisenforschung GmbH
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Dummy plasticity for purely elastic material
!--------------------------------------------------------------------------------------------------
submodule(phase:plastic) none

contains

!--------------------------------------------------------------------------------------------------
!> @brief Perform module initialization.
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function plastic_none_init() result(myPlasticity)

  logical, dimension(:), allocatable :: myPlasticity
  integer :: &
    ph
  type(tDict), pointer :: &
    phases


  myPlasticity = plastic_active('none')
  if (count(myPlasticity) == 0) return

  print'(/,1x,a)', '<<<+-  phase:mechanical:plastic:none init  -+>>>'
  print'(/,1x,a,1x,i0)', '# phases:',count(myPlasticity); flush(IO_STDOUT)


  phases => config_material%get_dict('phase')

  do ph = 1, phases%length
    if (.not. myPlasticity(ph)) cycle

    print'(/,1x,a,1x,i0,a)', 'phase',ph,': '//phases%key(ph)

    call phase_allocateState(plasticState(ph),count(material_ID_phase == ph),0,0,0)
  end do

end function plastic_none_init


end submodule none
