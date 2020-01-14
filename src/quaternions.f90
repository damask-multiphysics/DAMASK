!---------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @author Philip Eisenlohr, Michigan State University
!> @brief general quaternion math, not limited to unit quaternions
!> @details w is the real part, (x, y, z) are the imaginary parts.
!> @details https://en.wikipedia.org/wiki/Quaternion
!---------------------------------------------------------------------------------------------------
module quaternions
  use prec
  use IO

  implicit none
  public

  real(pReal), parameter, public :: P = -1.0_pReal                                                  !< parameter for orientation conversion.

  type, public :: quaternion
    real(pReal), private :: w = 0.0_pReal
    real(pReal), private :: x = 0.0_pReal
    real(pReal), private :: y = 0.0_pReal
    real(pReal), private :: z = 0.0_pReal


  contains
    procedure, private :: add__
    procedure, private :: pos__
    generic,   public  :: operator(+) => add__,pos__

    procedure, private :: sub__
    procedure, private :: neg__
    generic,   public  :: operator(-) => sub__,neg__

    procedure, private :: mul_quat__
    procedure, private :: mul_scal__
    generic,   public  :: operator(*) => mul_quat__, mul_scal__

    procedure, private :: div_quat__
    procedure, private :: div_scal__
    generic,   public  :: operator(/) => div_quat__, div_scal__

    procedure, private :: eq__
    generic,   public  :: operator(==) => eq__

    procedure, private :: neq__
    generic,   public  :: operator(/=) => neq__

    procedure, private :: pow_quat__
    procedure, private :: pow_scal__
    generic,   public  :: operator(**) => pow_quat__, pow_scal__

    procedure, public  :: abs   => abs__
    procedure, public  :: conjg => conjg__
    procedure, public  :: real  => real__
    procedure, public  :: aimag => aimag__

    procedure, public  :: homomorphed
    procedure, public  :: asArray
    procedure, public  :: inverse

  end type

  interface assignment (=)
    module procedure assign_quat__
    module procedure assign_vec__
  end interface assignment (=)
  
  interface quaternion
    module procedure init__
  end interface quaternion
  
  interface abs
    procedure abs__
  end interface abs
  
  interface dot_product
    procedure dot_product__
  end interface dot_product
  
  interface conjg
    module procedure conjg__
  end interface conjg
  
  interface exp
    module procedure exp__
  end interface exp
  
  interface log
    module procedure log__
  end interface log

  interface real
    module procedure real__
  end interface real

  interface aimag
    module procedure aimag__
  end interface aimag
  
  private :: &
    unitTest

contains


!--------------------------------------------------------------------------------------------------
!> @brief do self test
!--------------------------------------------------------------------------------------------------
subroutine quaternions_init

  write(6,'(/,a)') ' <<<+-  quaternions init  -+>>>'; flush(6)
  call unitTest

end subroutine quaternions_init


!---------------------------------------------------------------------------------------------------
!> construct a quaternion from a 4-vector
!---------------------------------------------------------------------------------------------------
type(quaternion) pure function init__(array)

  real(pReal), intent(in), dimension(4) :: array

  init__%w = array(1)
  init__%x = array(2)
  init__%y = array(3)
  init__%z = array(4)

end function init__


!---------------------------------------------------------------------------------------------------
!> assign a quaternion
!---------------------------------------------------------------------------------------------------
elemental pure subroutine assign_quat__(self,other)

  type(quaternion), intent(out) :: self
  type(quaternion), intent(in)  :: other

  self = [other%w,other%x,other%y,other%z]
  
end subroutine assign_quat__


!---------------------------------------------------------------------------------------------------
!> assign a 4-vector
!---------------------------------------------------------------------------------------------------
pure subroutine assign_vec__(self,other)

  type(quaternion), intent(out)                :: self
  real(pReal),       intent(in), dimension(4)  :: other

  self%w = other(1)
  self%x = other(2)
  self%y = other(3)
  self%z = other(4)

end subroutine assign_vec__


!---------------------------------------------------------------------------------------------------
!> add a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function add__(self,other)

  class(quaternion), intent(in) :: self,other

  add__ = [ self%w,  self%x,  self%y ,self%z] &
        + [other%w, other%x, other%y,other%z]
  
end function add__


!---------------------------------------------------------------------------------------------------
!> return (unary positive operator)
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function pos__(self)

  class(quaternion), intent(in) :: self

  pos__ = self * (+1.0_pReal)
  
end function pos__


!---------------------------------------------------------------------------------------------------
!> subtract a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function sub__(self,other)

  class(quaternion), intent(in) :: self,other

  sub__ = [ self%w,  self%x,  self%y ,self%z] &
        - [other%w, other%x, other%y,other%z]
  
end function sub__


!---------------------------------------------------------------------------------------------------
!> negate (unary negative operator)
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function neg__(self)

  class(quaternion), intent(in) :: self

  neg__ = self * (-1.0_pReal)
  
end function neg__


!---------------------------------------------------------------------------------------------------
!> multiply with a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function mul_quat__(self,other)

  class(quaternion), intent(in) :: self, other

  mul_quat__%w = self%w*other%w - self%x*other%x -      self%y*other%y - self%z*other%z
  mul_quat__%x = self%w*other%x + self%x*other%w + P * (self%y*other%z - self%z*other%y)
  mul_quat__%y = self%w*other%y + self%y*other%w + P * (self%z*other%x - self%x*other%z)
  mul_quat__%z = self%w*other%z + self%z*other%w + P * (self%x*other%y - self%y*other%x)

end function mul_quat__


!---------------------------------------------------------------------------------------------------
!> multiply with a scalar
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function mul_scal__(self,scal)

  class(quaternion), intent(in) :: self
  real(pReal),       intent(in) :: scal

  mul_scal__ = [self%w,self%x,self%y,self%z]*scal
  
end function mul_scal__


!---------------------------------------------------------------------------------------------------
!> divide by a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function div_quat__(self,other)

  class(quaternion), intent(in) :: self, other

  div_quat__ = self * (conjg(other)/(abs(other)**2.0_pReal))

end function div_quat__


!---------------------------------------------------------------------------------------------------
!> divide by a scalar
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function div_scal__(self,scal)

  class(quaternion), intent(in) :: self
  real(pReal),       intent(in) :: scal

  div_scal__ = [self%w,self%x,self%y,self%z]/scal

end function div_scal__


!---------------------------------------------------------------------------------------------------
!> test equality
!---------------------------------------------------------------------------------------------------
logical elemental pure function eq__(self,other)

  class(quaternion), intent(in) :: self,other

  eq__ = all(dEq([ self%w, self%x, self%y, self%z], &
                 [other%w,other%x,other%y,other%z]))

end function eq__


!---------------------------------------------------------------------------------------------------
!> test inequality
!---------------------------------------------------------------------------------------------------
logical elemental pure function neq__(self,other)

  class(quaternion), intent(in) :: self,other

  neq__ = .not. self%eq__(other)

end function neq__


!---------------------------------------------------------------------------------------------------
!> raise to the power of a quaternion
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function pow_quat__(self,expon)

  class(quaternion), intent(in) :: self
  type(quaternion),  intent(in) :: expon

  pow_quat__ = exp(log(self)*expon)

end function pow_quat__


!---------------------------------------------------------------------------------------------------
!> raise to the power of a scalar
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function pow_scal__(self,expon)

  class(quaternion), intent(in) :: self
  real(pReal),       intent(in) :: expon

  pow_scal__ = exp(log(self)*expon)

end function pow_scal__


!---------------------------------------------------------------------------------------------------
!> take exponential
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function exp__(a)

  class(quaternion), intent(in) :: a
  real(pReal)                   :: absImag

  absImag = norm2(aimag(a))

  exp__ = merge(exp(a%w) * [               cos(absImag), &
                             a%x/absImag * sin(absImag), &
                             a%y/absImag * sin(absImag), &
                             a%z/absImag * sin(absImag)], &
                IEEE_value(1.0_pReal,IEEE_SIGNALING_NAN), &
                dNeq0(absImag))

end function exp__


!---------------------------------------------------------------------------------------------------
!> take logarithm
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function log__(a)

  class(quaternion), intent(in) :: a
  real(pReal)                   :: absImag

  absImag = norm2(aimag(a))

  log__ = merge([log(abs(a)), &
                 a%x/absImag * acos(a%w/abs(a)), &
                 a%y/absImag * acos(a%w/abs(a)), &
                 a%z/absImag * acos(a%w/abs(a))], &
                IEEE_value(1.0_pReal,IEEE_SIGNALING_NAN), &
                dNeq0(absImag))

end function log__


!---------------------------------------------------------------------------------------------------
!> return norm
!---------------------------------------------------------------------------------------------------
real(pReal) elemental pure function abs__(self)

  class(quaternion), intent(in) :: self

  abs__ = norm2([self%w,self%x,self%y,self%z])

end function abs__


!---------------------------------------------------------------------------------------------------
!> calculate dot product
!---------------------------------------------------------------------------------------------------
real(pReal) elemental pure function dot_product__(a,b)

  class(quaternion), intent(in) :: a,b

  dot_product__ = a%w*b%w + a%x*b%x + a%y*b%y + a%z*b%z

end function dot_product__


!---------------------------------------------------------------------------------------------------
!> take conjugate complex
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function conjg__(self)

  class(quaternion), intent(in) :: self

  conjg__ = [self%w,-self%x,-self%y,-self%z]

end function conjg__


!---------------------------------------------------------------------------------------------------
!> homomorph
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function homomorphed(self)

  class(quaternion), intent(in) :: self

  homomorphed = - self

end function homomorphed


!---------------------------------------------------------------------------------------------------
!> return as plain array
!---------------------------------------------------------------------------------------------------
pure function asArray(self)

  real(pReal), dimension(4)     :: asArray
  class(quaternion), intent(in) :: self

  asArray = [self%w,self%x,self%y,self%z]

end function asArray


!---------------------------------------------------------------------------------------------------
!> real part (scalar)
!---------------------------------------------------------------------------------------------------
pure function real__(self)

  real(pReal)                   :: real__
  class(quaternion), intent(in) :: self

  real__ = self%w

end function real__


!---------------------------------------------------------------------------------------------------
!> imaginary part (3-vector)
!---------------------------------------------------------------------------------------------------
pure function aimag__(self)

  real(pReal), dimension(3)     :: aimag__
  class(quaternion), intent(in) :: self

  aimag__ = [self%x,self%y,self%z]

end function aimag__


!---------------------------------------------------------------------------------------------------
!> inverse
!---------------------------------------------------------------------------------------------------
type(quaternion) elemental pure function inverse(self)

  class(quaternion), intent(in) :: self

  inverse = conjg(self)/abs(self)**2.0_pReal

end function inverse


!--------------------------------------------------------------------------------------------------
!> @brief check correctness of (some) quaternions functions
!--------------------------------------------------------------------------------------------------
subroutine unitTest

  real(pReal), dimension(4) :: qu
  type(quaternion)          :: q, q_2

  call random_number(qu)
  qu = (qu-0.5_pReal) * 2.0_pReal
  q  = quaternion(qu)

  q_2= qu
  if(any(dNeq(q%asArray(),q_2%asArray())))             call IO_error(401,ext_msg='assign_vec__')

  q_2 = q + q
  if(any(dNeq(q_2%asArray(),2.0_pReal*qu)))            call IO_error(401,ext_msg='add__')

  q_2 = q - q
  if(any(dNeq0(q_2%asArray())))                        call IO_error(401,ext_msg='sub__')

  q_2 = q * 5.0_pReal
  if(any(dNeq(q_2%asArray(),5.0_pReal*qu)))            call IO_error(401,ext_msg='mul__')

  q_2 = q / 0.5_pReal
  if(any(dNeq(q_2%asArray(),2.0_pReal*qu)))            call IO_error(401,ext_msg='div__')

  q_2 = q * 0.3_pReal
  if(dNeq0(abs(q)) .and. q_2 == q)                     call IO_error(401,ext_msg='eq__')

  q_2 = q
  if(q_2 /= q)                                         call IO_error(401,ext_msg='neq__')

  if(dNeq(abs(q),norm2(qu)))                           call IO_error(401,ext_msg='abs__')
  if(dNeq(abs(q)**2.0_pReal, real(q*q%conjg()),1.0e-14_pReal)) &
                                                       call IO_error(401,ext_msg='abs__/*conjg')

  if(any(dNeq(q%asArray(),qu)))                        call IO_error(401,ext_msg='eq__')
  if(dNeq(q%real(),       qu(1)))                      call IO_error(401,ext_msg='real()')
  if(any(dNeq(q%aimag(),  qu(2:4))))                   call IO_error(401,ext_msg='aimag()')

  q_2 = q%homomorphed()
  if(q                 /= q_2*    (-1.0_pReal))        call IO_error(401,ext_msg='homomorphed')
  if(dNeq(q_2%real(),     qu(1)*  (-1.0_pReal)))       call IO_error(401,ext_msg='homomorphed/real')
  if(any(dNeq(q_2%aimag(),qu(2:4)*(-1.0_pReal))))      call IO_error(401,ext_msg='homomorphed/aimag')

  q_2 = conjg(q)
  if(dNeq(abs(q),abs(q_2)))                            call IO_error(401,ext_msg='conjg/abs')
  if(q /= conjg(q_2))                                  call IO_error(401,ext_msg='conjg/involution')
  if(dNeq(q_2%real(),     q%real()))                   call IO_error(401,ext_msg='conjg/real')
  if(any(dNeq(q_2%aimag(),q%aimag()*(-1.0_pReal))))    call IO_error(401,ext_msg='conjg/aimag')

  if(abs(q) > 0.0_pReal) then
    q_2 = q * q%inverse()
    if(     dNeq(real(q_2), 1.0_pReal,1.0e-15_pReal))  call IO_error(401,ext_msg='inverse/real')
    if(any(dNeq0(aimag(q_2),          1.0e-15_pReal))) call IO_error(401,ext_msg='inverse/aimag')

    q_2 = q/abs(q)
    q_2 = conjg(q_2) - inverse(q_2)
    if(any(dNeq0(q_2%asArray(),1.0e-15_pReal)))        call IO_error(401,ext_msg='inverse/conjg')
  endif

#if !(defined(__GFORTRAN__) &&  __GNUC__ < 9)
  if (norm2(aimag(q)) > 0.0_pReal) then
    if (dNeq0(abs(q-exp(log(q))),1.0e-13_pReal))       call IO_error(401,ext_msg='exp/log')
    if (dNeq0(abs(q-log(exp(q))),1.0e-13_pReal))       call IO_error(401,ext_msg='log/exp')
  endif
#endif

end subroutine unitTest


end module quaternions
