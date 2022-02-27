!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Dummy homogenization scheme for 1 constituent per material point
!--------------------------------------------------------------------------------------------------
submodule(homogenization:damage) damage_pass

contains

module subroutine pass_init()

  integer :: &
    ho

  print'(/,1x,a)', '<<<+-  homogenization:damage:pass init  -+>>>'

  do ho = 1, size(damage_active)

    if (.not. damage_active(ho)) cycle

    if (homogenization_Nconstituents(ho) /= 1) &
      call IO_error(211,ext_msg='(pass) with N_constituents !=1')
  end do

end subroutine pass_init

end submodule damage_pass
