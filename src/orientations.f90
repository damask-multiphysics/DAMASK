module orientations
 use rotations
 use prec, only: &
   pStringLen

 implicit none
 type, extends(rotation), public :: orientation
   character(len=pStringLen) :: sym = 'none'
 end type orientation

 interface orientation
   module procedure :: orientation_init
 end interface orientation

contains

type(orientation) function orientation_init(sym,eu,ax,om,qu,cu,ho,ro)
  use prec
  implicit none
  character(len=pStringLen), intent(in), optional                 :: sym
  real(pReal),               intent(in), optional, dimension(3)   :: eu, cu, ho
  real(pReal),               intent(in), optional, dimension(4)   :: ax, qu, ro
  real(pReal),               intent(in), optional, dimension(3,3) :: om

  if (present(sym)) orientation_init%sym = sym

  if (present(om)) then
    call orientation_init%fromRotationMatrix(om)
  endif

end function orientation_init

end module
