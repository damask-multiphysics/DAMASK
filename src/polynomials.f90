!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Polynomial representation for variable data.
!--------------------------------------------------------------------------------------------------
module polynomials
  use prec
  use IO
  use YAML_parse
  use YAML_types

  implicit none(type,external)
  private

  type, public :: tPolynomial
    real(pREAL), dimension(:), allocatable :: coef
    real(pREAL) :: x_ref = huge(0.0_pREAL)
    contains
    procedure, public :: at => eval
  end type tPolynomial

  interface polynomial
    module procedure polynomial_from_coef
    module procedure polynomial_from_dict
  end interface polynomial

  public :: &
    polynomial, &
    polynomials_init

contains


!--------------------------------------------------------------------------------------------------
!> @brief Run self-test.
!--------------------------------------------------------------------------------------------------
subroutine polynomials_init()

  print'(/,1x,a)', '<<<+-  polynomials init  -+>>>'; flush(IO_STDOUT)

  call selfTest()

end subroutine polynomials_init


!--------------------------------------------------------------------------------------------------
!> @brief Initialize a polynomial from coefficients.
!--------------------------------------------------------------------------------------------------
pure function polynomial_from_coef(coef,x_ref) result(p)

  real(pREAL), dimension(0:), intent(in) :: coef
  real(pREAL), intent(in) :: x_ref
  type(tPolynomial) :: p


  p%coef = coef
  p%x_ref = x_ref

end function polynomial_from_coef


!--------------------------------------------------------------------------------------------------
!> @brief Initialize a polynomial from a dictionary with coefficients.
!--------------------------------------------------------------------------------------------------
function polynomial_from_dict(dict,y,x) result(p)

  type(tDict), intent(in) :: dict
  character(len=*), intent(in) :: y, x
  type(tPolynomial) :: p

  real(pREAL), dimension(:), allocatable :: coef
  real(pREAL) :: x_ref
  integer :: i, o
  character(len=1) :: o_s


  allocate(coef(1),source=dict%get_asReal(y))

  if (dict%contains(y//','//x)) then
    x_ref = dict%get_asReal(x//'_ref')
    coef = [coef,dict%get_asReal(y//','//x)]
  end if
  do o = 2,4
    write(o_s,'(I0.0)') o
    if (dict%contains(y//','//x//'^'//o_s)) then
      x_ref = dict%get_asReal(x//'_ref')
      coef = [coef,[(0.0_pREAL,i=size(coef),o-1)],dict%get_asReal(y//','//x//'^'//o_s)]
    end if
  end do

  p = Polynomial(coef,x_ref)

end function polynomial_from_dict


!--------------------------------------------------------------------------------------------------
!> @brief Evaluate a polynomial.
!> @details https://nvlpubs.nist.gov/nistpubs/jres/71b/jresv71bn1p11_a1b.pdf (eq. 1.2)
!--------------------------------------------------------------------------------------------------
pure function eval(self,x) result(y)

  class(tPolynomial), intent(in) :: self
  real(pREAL), intent(in) :: x
  real(pREAL) :: y

  integer :: o


  y = self%coef(ubound(self%coef,1))
  do o = ubound(self%coef,1)-1, 0, -1
#ifndef __INTEL_LLVM_COMPILER
    y = y*(x-self%x_ref) +self%coef(o)
#else
    y = IEEE_FMA(y,x-self%x_ref,self%coef(o))
#endif
  end do

end function eval


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of polynomical functionality.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  type(tPolynomial) :: p1, p2
  real(pREAL), dimension(5) :: coef
  integer :: i
  real(pREAL) :: x_ref, x, y
  type(tDict), pointer :: dict
  character(len=pSTRLEN), dimension(size(coef)) :: coef_s
  character(len=pSTRLEN) :: x_ref_s, x_s, YAML_s


  call random_number(coef)
  call random_number(x_ref)
  call random_number(x)

  coef = coef*10_pREAL -0.5_pREAL
  x_ref = x_ref*10_pREAL -0.5_pREAL
  x = x*10_pREAL -0.5_pREAL

  p1 = polynomial([coef(1)],x_ref)
  if (dNeq(p1%at(x),coef(1)))      error stop 'polynomial: eval(constant)'

  p1 = polynomial(coef,x_ref)
  if (dNeq(p1%at(x_ref),coef(1)))  error stop 'polynomial: @ref'

  do i = 1, size(coef_s)
    write(coef_s(i),*) coef(i)
  end do
  write(x_ref_s,*) x_ref
  write(x_s,*) x
  YAML_s = 'C: '//trim(adjustl(coef_s(1)))//IO_EOL//&
           'C,T: '//trim(adjustl(coef_s(2)))//IO_EOL//&
           'C,T^2: '//trim(adjustl(coef_s(3)))//IO_EOL//&
           'C,T^3: '//trim(adjustl(coef_s(4)))//IO_EOL//&
           'C,T^4: '//trim(adjustl(coef_s(5)))//IO_EOL//&
           'T_ref: '//trim(adjustl(x_ref_s))//IO_EOL
  dict => YAML_parse_str_asDict(trim(YAML_s))
  p2 = polynomial(dict,'C','T')
  if (dNeq(p1%at(x),p2%at(x),1.0e-6_pREAL))                      error stop 'polynomials: init'
  y = coef(1)+coef(2)*(x-x_ref)+coef(3)*(x-x_ref)**2+coef(4)*(x-x_ref)**3+coef(5)*(x-x_ref)**4
  if (dNeq(p1%at(x),y,1.0e-6_pREAL))                             error stop 'polynomials: eval(full)'

  YAML_s = 'C: 0.0'//IO_EOL//&
           'C,T: '//trim(adjustl(coef_s(2)))//IO_EOL//&
           'T_ref: '//trim(adjustl(x_ref_s))//IO_EOL
  dict => YAML_parse_str_asDict(trim(YAML_s))
  p1 = polynomial(dict,'C','T')
  if (dNeq(p1%at(x_ref+x),-p1%at(x_ref-x),1.0e-10_pREAL))         error stop 'polynomials: eval(linear)'

  YAML_s = 'C: 0.0'//IO_EOL//&
           'C,T^2: '//trim(adjustl(coef_s(3)))//IO_EOL//&
           'T_ref: '//trim(adjustl(x_ref_s))//IO_EOL
  dict => YAML_parse_str_asDict(trim(YAML_s))
  p1 = polynomial(dict,'C','T')
  if (dNeq(p1%at(x_ref+x),p1%at(x_ref-x),1e-10_pREAL))            error stop 'polynomials: eval(quadratic)'

  YAML_s = 'Y: '//trim(adjustl(coef_s(1)))//IO_EOL//&
           'Y,X^3: '//trim(adjustl(coef_s(2)))//IO_EOL//&
           'X_ref: '//trim(adjustl(x_ref_s))//IO_EOL
  dict => YAML_parse_str_asDict(trim(YAML_s))
  p1 = polynomial(dict,'Y','X')
  if (dNeq(p1%at(x_ref+x)-coef(1),-(p1%at(x_ref-x)-coef(1)),1.0e-8_pREAL)) error stop 'polynomials: eval(cubic)'

  YAML_s = 'Y: '//trim(adjustl(coef_s(1)))//IO_EOL//&
           'Y,X^4: '//trim(adjustl(coef_s(2)))//IO_EOL//&
           'X_ref: '//trim(adjustl(x_ref_s))//IO_EOL
  dict => YAML_parse_str_asDict(trim(YAML_s))
  p1 = polynomial(dict,'Y','X')
  if (dNeq(p1%at(x_ref+x),p1%at(x_ref-x),1.0e-6_pREAL))           error stop 'polynomials: eval(quartic)'


end subroutine selfTest

end module polynomials
