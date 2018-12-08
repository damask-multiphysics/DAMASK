module orientations
 use rotations

 implicit none
 type, extends(rotation), public :: orientation
 end type orientation

 interface orientation
   module procedure :: orientation_init
 end interface orientation

contains

type(orientation) function orientation_init(eu,ax,om,qu,cu,ho,ro)
  use prec
  implicit none
  real(pReal),      intent(in), optional, dimension(3)   :: eu, cu, ho
  real(pReal),      intent(in), optional, dimension(4)   :: ax, qu, ro
  real(pReal),      intent(in), optional, dimension(3,3) :: om
 
  if (present(om)) then
    call orientation_init%fromRotationMatrix(om)
  endif

end function orientation_init

end module
