!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, KU Leuven
!> @brief Polynomial representation for variable data
!--------------------------------------------------------------------------------------------------
module polynomials
  use prec
  use IO
  use YAML_parse
  use YAML_types

  implicit none
  private

  type, public :: tPolynomial
    real(pReal), dimension(:), allocatable :: coef
    real(pReal) :: x_ref
    contains
    procedure, public :: at => eval
    procedure, public :: der1_at => eval_der1
  end type tPolynomial

  interface polynomial
    module procedure polynomial_from_dict
    module procedure polynomial_from_coef
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
!> @brief Initialize a Polynomial from Coefficients.
!--------------------------------------------------------------------------------------------------
function polynomial_from_coef(coef,x_ref) result(p)

  real(pReal), dimension(:), intent(in) :: coef
  real(pReal), intent(in) :: x_ref
  type(tPolynomial) :: p


  allocate(p%coef(0:size(coef)-1),source=coef)                                                      ! should be zero based
  p%x_ref = x_ref

end function polynomial_from_coef


!--------------------------------------------------------------------------------------------------
!> @brief Initialize a Polynomial from a Dictionary with Coefficients.
!--------------------------------------------------------------------------------------------------
function polynomial_from_dict(dict,y,x) result(p)

  type(tDict), intent(in) :: dict
  character(len=*), intent(in) :: y, x
  type(tPolynomial) :: p

  real(pReal), dimension(:), allocatable :: coef
  real(pReal) :: x_ref


  allocate(coef(1),source=dict%get_asFloat(y))

  if (dict%contains(y//','//x)) then
    x_ref = dict%get_asFloat(x//'_ref')
    coef = [coef,dict%get_asFloat(y//','//x)]
    if (dict%contains(y//','//x//'^2')) then
      coef = [coef,dict%get_asFloat(y//','//x//'^2')]
    end if
  else
    x_ref = huge(0.0_pReal)                                                                         ! Simplify debugging
  end if

  p = Polynomial(coef,x_ref)

end function polynomial_from_dict


!--------------------------------------------------------------------------------------------------
!> @brief Evaluate a Polynomial.
!--------------------------------------------------------------------------------------------------
pure function eval(self,x) result(y)

  class(tPolynomial), intent(in) :: self
  real(pReal), intent(in) :: x
  real(pReal) :: y

  integer :: i


  y = self%coef(0)
  do i = 1, ubound(self%coef,1)
    y = y + self%coef(i) * (x-self%x_ref)**i
  enddo

end function eval


!--------------------------------------------------------------------------------------------------
!> @brief Evaluate a first derivative of Polynomial.
!--------------------------------------------------------------------------------------------------
pure function eval_der1(self,x) result(y)

  class(tPolynomial), intent(in) :: self
  real(pReal), intent(in) :: x
  real(pReal) :: y

  integer :: i


  y = 0.0_pReal
  do i = 1, ubound(self%coef,1)
    y = y + real(i,pReal)*self%coef(i) * (x-self%x_ref)**(i-1)
  enddo

end function eval_der1


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of polynomical functionality.
!--------------------------------------------------------------------------------------------------
subroutine selfTest

  type(tPolynomial) :: p1, p2
  real(pReal), dimension(3) :: coef
  real(pReal) :: x_ref, x
  class(tNode), pointer :: dict
  character(len=pStringLen), dimension(3) :: coef_s
  character(len=pStringLen) :: x_ref_s, x_s, YAML_s

  call random_number(coef)
  call random_number(x_ref)
  call random_number(x)

  coef = coef*10_pReal -0.5_pReal
  x_ref = x_ref*10_pReal -0.5_pReal
  x = x*10_pReal -0.5_pReal

  p1 = polynomial(coef,x_ref)
  if (dNeq(p1%at(x_ref),coef(1)))  error stop 'polynomial: @ref'

  write(coef_s(1),*) coef(1)
  write(coef_s(2),*) coef(2)
  write(coef_s(3),*) coef(3)
  write(x_ref_s,*) x_ref
  write(x_s,*) x
  YAML_s = 'C: '//trim(adjustl(coef_s(1)))//IO_EOL//&
           'C,T: '//trim(adjustl(coef_s(2)))//IO_EOL//&
           'C,T^2: '//trim(adjustl(coef_s(3)))//IO_EOL//&
           'T_ref: '//trim(adjustl(x_ref_s))//IO_EOL
  Dict => YAML_parse_str(trim(YAML_s))
  p2 = polynomial(dict%asDict(),'C','T')
  if (dNeq(p1%at(x),p2%at(x),1.0e-10_pReal))                      error stop 'polynomials: init'

  p1 = polynomial(coef*[0.0_pReal,1.0_pReal,0.0_pReal],x_ref)
  if (dNeq(p1%at(x_ref+x),-p1%at(x_ref-x),1.0e-10_pReal))         error stop 'polynomials: eval(odd)'
  if (dNeq(p1%der1_at(x),p1%der1_at(5.0_pReal*x),1.0e-10_pReal))  error stop 'polynomials: eval_der(odd)'

  p1 = polynomial(coef*[0.0_pReal,0.0_pReal,1.0_pReal],x_ref)
  if (dNeq(p1%at(x_ref+x),p1%at(x_ref-x),1e-10_pReal))            error stop 'polynomials: eval(even)'
  if (dNeq(p1%der1_at(x_ref+x),-p1%der1_at(x_ref-x),1e-10_pReal)) error stop 'polynomials: eval_der(even)'


end subroutine selfTest

end module polynomials
