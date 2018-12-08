! ###################################################################
! Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
!     - Redistributions of source code must retain the above copyright notice, this list 
!        of conditions and the following disclaimer.
!     - Redistributions in binary form must reproduce the above copyright notice, this 
!        list of conditions and the following disclaimer in the documentation and/or 
!        other materials provided with the distribution.
!     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
!        of its contributors may be used to endorse or promote products derived from 
!        this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
! CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
! OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
! USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
! ###################################################################

module quaternions
 
 use prec
 implicit none

 public
 type, public :: quaternion
   real(pReal) :: w = 0.0_pReal
   real(pReal) :: x = 0.0_pReal
   real(pReal) :: y = 0.0_pReal
   real(pReal) :: z = 0.0_pReal
 

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

   procedure, private :: abs__
   procedure, private :: dot_product__
   procedure, private :: conjg__
   procedure, private :: exp__
   procedure, private :: log__

   procedure, public  :: homomorphed => quat_homomorphed

   !procedure,private  :: quat_write
   !generic :: write(formatted) => quat_write

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

contains


!--------------------------------------------------------------------------
!> constructor for a quaternion from a 4-vector
!--------------------------------------------------------------------------
type(quaternion) pure function init__(array)

 implicit none
 real(pReal), intent(in), dimension(4) :: array
    
 init__%w=array(1)
 init__%x=array(2)
 init__%y=array(3)
 init__%z=array(4)

end function init__


!--------------------------------------------------------------------------
!> assing a quaternion
!--------------------------------------------------------------------------
elemental subroutine assign_quat__(self,other)

 implicit none
 type(quaternion), intent(out) :: self
 type(quaternion), intent(in)  :: other
    
 self%w = other%w
 self%x = other%x
 self%y = other%y
 self%z = other%z
  
end subroutine assign_quat__


!--------------------------------------------------------------------------
!> assing a 4-vector
!--------------------------------------------------------------------------
pure subroutine assign_vec__(self,other)

 implicit none
 type(quaternion), intent(out)                :: self
 real(pReal),         intent(in), dimension(4)  :: other
    
 self%w = other(1)
 self%x = other(2)
 self%y = other(3)
 self%z = other(4)
  
end subroutine assign_vec__


!--------------------------------------------------------------------------
!> addition of two quaternions
!--------------------------------------------------------------------------
type(quaternion) elemental function add__(self,other)

 implicit none
 class(quaternion), intent(in) :: self,other
    
 add__%w = self%w + other%w
 add__%x = self%x + other%x
 add__%y = self%y + other%y
 add__%z = self%z + other%z
  
end function add__


!--------------------------------------------------------------------------
!> unary positive operator
!--------------------------------------------------------------------------
type(quaternion) elemental function pos__(self)

 implicit none
 class(quaternion), intent(in) :: self
    
 pos__%w = self%w 
 pos__%x = self%x 
 pos__%y = self%y 
 pos__%z = self%z 
  
end function pos__


!--------------------------------------------------------------------------
!> subtraction of two quaternions
!--------------------------------------------------------------------------
type(quaternion) elemental function sub__(self,other)

 implicit none
 class(quaternion), intent(in) :: self,other
    
 sub__%w = self%w - other%w
 sub__%x = self%x - other%x
 sub__%y = self%y - other%y
 sub__%z = self%z - other%z
  
end function sub__


!--------------------------------------------------------------------------
!> unary positive operator
!--------------------------------------------------------------------------
type(quaternion) elemental function neg__(self)

 implicit none
 class(quaternion), intent(in) :: self
    
 neg__%w = -self%w 
 neg__%x = -self%x 
 neg__%y = -self%y 
 neg__%z = -self%z 
  
end function neg__


!--------------------------------------------------------------------------
!> multiplication of two quaternions
!--------------------------------------------------------------------------
type(quaternion) elemental function mul_quat__(self,other)

 implicit none
 class(quaternion), intent(in) :: self, other

 mul_quat__%w = self%w*other%w - self%x*other%x -           self%y*other%y - self%z*other%z
 mul_quat__%x = self%w*other%x + self%x*other%w + epsijk * (self%y*other%z - self%z*other%y)
 mul_quat__%y = self%w*other%y + self%y*other%w + epsijk * (self%z*other%x - self%x*other%z)
 mul_quat__%z = self%w*other%z + self%z*other%w + epsijk * (self%x*other%y - self%y*other%x)
    
end function mul_quat__


!--------------------------------------------------------------------------
!> multiplication of quaternions with scalar
!--------------------------------------------------------------------------
type(quaternion) elemental function mul_scal__(self,scal)

 implicit none
 class(quaternion), intent(in) :: self
 real(pReal), intent(in) :: scal

 mul_scal__%w = self%w*scal
 mul_scal__%x = self%x*scal
 mul_scal__%y = self%y*scal
 mul_scal__%z = self%z*scal
    
end function mul_scal__


!--------------------------------------------------------------------------
!> division of two quaternions
!--------------------------------------------------------------------------
type(quaternion) elemental function div_quat__(self,other)

 implicit none
 class(quaternion), intent(in) :: self, other

 div_quat__ = self * (conjg(other)/(abs(other)**2.0_pReal))

end function div_quat__


!--------------------------------------------------------------------------
!> divisiont of quaternions by scalar
!--------------------------------------------------------------------------
type(quaternion) elemental function div_scal__(self,scal)

 implicit none
 class(quaternion), intent(in) :: self
 real(pReal), intent(in) :: scal

 div_scal__ = [self%w,self%x,self%y,self%z]/scal

end function div_scal__


!--------------------------------------------------------------------------
!> equality of two quaternions
!--------------------------------------------------------------------------
logical elemental function eq__(self,other)
 implicit none
 class(quaternion), intent(in) :: self,other

 eq__ = all(dEq([ self%w, self%x, self%y, self%z], &
                [other%w,other%x,other%y,other%z]))
    
end function eq__


!--------------------------------------------------------------------------
!> inequality of two quaternions
!--------------------------------------------------------------------------
logical elemental function neq__(self,other)

 implicit none
 class(quaternion), intent(in) :: self,other

 neq__ = .not. self%eq__(other)
    
end function neq__


!--------------------------------------------------------------------------
!> quaternion to the power of a scalar
!--------------------------------------------------------------------------
type(quaternion) elemental function pow_scal__(self,expon)

 implicit none
 class(quaternion), intent(in) :: self
 real(pReal), intent(in) :: expon
 
 pow_scal__ = exp(log(self)*expon)
    
end function pow_scal__


!--------------------------------------------------------------------------
!> quaternion to the power of a quaternion
!--------------------------------------------------------------------------
type(quaternion) elemental function pow_quat__(self,expon)

 implicit none
 class(quaternion), intent(in) :: self
 type(quaternion),  intent(in) :: expon
 
 pow_quat__ = exp(log(self)*expon)
    
end function pow_quat__


!--------------------------------------------------------------------------
!> exponential of a quaternion
!> ToDo: Lacks any check for invalid operations
!--------------------------------------------------------------------------
type(quaternion) elemental function exp__(self)

 implicit none
 class(quaternion), intent(in) :: self
 real(pReal)                     :: absImag

 absImag = norm2([self%x, self%y, self%z])

 exp__ = exp(self%w) * [                 cos(absImag), &
                        self%x/absImag * sin(absImag), &
                        self%y/absImag * sin(absImag), &
                        self%z/absImag * sin(absImag)]

end function exp__


!--------------------------------------------------------------------------
!> logarithm of a quaternion
!> ToDo: Lacks any check for invalid operations
!--------------------------------------------------------------------------
type(quaternion) elemental function log__(self)

 implicit none
 class(quaternion), intent(in) :: self
 real(pReal)                     :: absImag

 absImag = norm2([self%x, self%y, self%z])

 log__ = [log(abs(self)), &
          self%x/absImag * acos(self%w/abs(self)), &
          self%y/absImag * acos(self%w/abs(self)), &
          self%z/absImag * acos(self%w/abs(self))]
    
end function log__


!--------------------------------------------------------------------------
!> norm of a quaternion
!--------------------------------------------------------------------------
real(pReal) elemental function abs__(a)

 implicit none
 class(quaternion), intent(in) :: a

 abs__ = norm2([a%w,a%x,a%y,a%z])
    
end function abs__


!--------------------------------------------------------------------------
!> dot product of two quaternions
!--------------------------------------------------------------------------
real(pReal) elemental function dot_product__(a,b)

 implicit none
 class(quaternion), intent(in) :: a,b

 dot_product__ = a%w*b%w + a%x*b%x + a%y*b%y + a%z*b%z
    
end function dot_product__


!--------------------------------------------------------------------------
!> conjugate complex of a quaternion
!--------------------------------------------------------------------------
type(quaternion) elemental function conjg__(a)

 implicit none
 class(quaternion), intent(in) :: a

 conjg__ = quaternion([a%w, -a%x, -a%y, -a%z])
    
end function conjg__


!--------------------------------------------------------------------------
!> homomorphed quaternion of a quaternion
!--------------------------------------------------------------------------
type(quaternion) elemental function quat_homomorphed(a)

 implicit none
 class(quaternion), intent(in) :: a

 quat_homomorphed = quaternion(-[a%w,a%x,a%y,a%z])
    
end function quat_homomorphed

end module quaternions
