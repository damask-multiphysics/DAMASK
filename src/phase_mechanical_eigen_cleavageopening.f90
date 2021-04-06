!--------------------------------------------------------------------------------------------------
!> @author Luv Sharma, Max-Planck-Institut für Eisenforschung GmbH
!> @author Pratheek Shanthraj, Max-Planck-Institut für Eisenforschung GmbH
!> @brief material subroutine incorporating kinematics resulting from opening of cleavage planes
!> @details to be done
!--------------------------------------------------------------------------------------------------
submodule(phase:eigen) cleavageopening

contains


!--------------------------------------------------------------------------------------------------
!> @brief module initialization
!> @details reads in material parameters, allocates arrays, and does sanity checks
!--------------------------------------------------------------------------------------------------
module function damage_anisobrittle_init() result(myKinematics)

  logical, dimension(:), allocatable :: myKinematics


  myKinematics = kinematics_active2('anisobrittle')
  if(count(myKinematics) == 0) return

  print'(/,a)', ' <<<+-  phase:mechanical:eigen:cleavageopening init  -+>>>'
  print'(a,i2)', ' # phases: ',count(myKinematics); flush(IO_STDOUT)

end function damage_anisobrittle_init


end submodule cleavageopening
