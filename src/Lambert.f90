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

!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!
!> @brief everything that has to do with the modified Lambert projections
!
!> @details This module contains a number of projection functions for the modified 
!> Lambert projection between square lattice and 2D hemisphere, hexagonal lattice
!> and 2D hemisphere, as well as the more complex mapping between a 3D cubic grid
!> and the unit quaternion hemisphere with positive scalar component.  In addition, there 
!> are some other projections, such as the stereographic one.  Each function is named
!> by the projection, the dimensionality of the starting grid, and the forward or inverse
!> character.  For each function, there is also a single precision and a double precision
!> version, but we use the interface formalism to have only a single call.  The Forward
!> mapping is taken to be the one from the simple grid to the curved grid.  Since the module
!> deals with various grids, we also add a few functions/subroutines that apply symmetry
!> operations on those grids.
!> References:
!> D. Rosca, A. Morawiec, and M. De Graef. “A new method of constructing a grid 
!> in the space of 3D rotations and its applications to texture analysis”. 
!> Modeling and Simulations in Materials Science and Engineering 22, 075013 (2014).
!--------------------------------------------------------------------------
module Lambert
 use math
 use prec

 implicit none

 real(pReal), private :: &
   sPi  = sqrt(PI), &
   pref = sqrt(6.0_pReal/PI), &
   ! the following constants are used for the cube to quaternion hemisphere mapping
   ap   = PI**(2.0_pReal/3.0_pReal), &
   sc   = 0.897772786961286_pReal, &      ! a/ap
   beta = 0.962874509979126_pReal, &      ! pi^(5/6)/6^(1/6)/2
   R1   = 1.330670039491469_pReal, &      ! (3pi/4)^(1/3)
   r2   = sqrt(2.0_pReal), &
   pi12 = PI/12.0_pReal, &
   prek = 1.643456402972504_pReal, &      ! R1 2^(1/4)/beta
   r24  = sqrt(24.0_pReal)

 private
 public :: &
   LambertCubeToBall, &
   LambertBallToCube
 private :: &
  GetPyramidOrder

contains


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief map from 3D cubic grid to 3D ball 
!--------------------------------------------------------------------------
function LambertCubeToBall(cube) result(ball)
 use, intrinsic :: IEEE_ARITHMETIC

 implicit none
 real(pReal), intent(in), dimension(3) :: cube
 real(pReal),             dimension(3) :: ball, LamXYZ, XYZ
 real(pReal)                  ::  T(2), c, s, q
 real(pReal), parameter       :: eps = 1.0e-8_pReal
 integer(pInt), dimension(3) :: p
 integer(pInt), dimension(2) :: order
 
 if (maxval(abs(cube)) > ap/2.0+eps) then
   ball = IEEE_value(cube,IEEE_positive_inf)
   return
 end if
 
 ! transform to the sphere grid via the curved square, and intercept the zero point
 center: if (all(dEq0(cube))) then
   ball  = 0.0_pReal
 else center
   ! get pyramide and scale by grid parameter ratio
   p = GetPyramidOrder(cube)
   XYZ = cube(p) * sc
 
   ! intercept all the points along the z-axis
   special: if (all(dEq0(XYZ(1:2)))) then
     LamXYZ = [ 0.0_pReal, 0.0_pReal, pref * XYZ(3) ]
   else special
     order = merge( [2,1], [1,2], abs(XYZ(2)) <= abs(XYZ(1)))                                        ! order of absolute values of XYZ
     q = pi12 *  XYZ(order(1))/XYZ(order(2))                                                         ! smaller by larger
     c = cos(q)
     s = sin(q)
     q = prek * XYZ(order(2))/ sqrt(r2-c)
     T = [ (r2*c - 1.0), r2 * s] * q
 
 ! transform to sphere grid (inverse Lambert)
 ! [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
     c = sum(T**2)
     s = Pi * c/(24.0*XYZ(3)**2)
     c = sPi * c / r24 / XYZ(3)
     q = sqrt( 1.0 - s )
     LamXYZ = [ T(order(2)) * q, T(order(1)) * q, pref * XYZ(3) - c ]
   endif special
 
 ! reverse the coordinates back to the regular order according to the original pyramid number
   ball = LamXYZ(p)

 endif center

end function LambertCubeToBall

!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief map from 3D ball to 3D cubic grid  
!--------------------------------------------------------------------------
pure function LambertBallToCube(xyz) result(cube)
 use, intrinsic :: IEEE_ARITHMETIC

 implicit none
 real(pReal), intent(in), dimension(3) :: xyz
 real(pReal),             dimension(3) :: cube, xyz1, xyz3
 real(pReal),             dimension(2) :: Tinv, xyz2
 real(pReal)                           :: rs, qxy, q2, sq2, q, tt
 integer(pInt)         , dimension(3) :: p
 
 rs = norm2(xyz)
 if (rs > R1) then 
   cube = IEEE_value(cube,IEEE_positive_inf)
   return
 endif
 
 center: if (all(dEq0(xyz))) then 
   cube = 0.0_pReal
 else center
   p = GetPyramidOrder(xyz)
   xyz3 = xyz(p)
   
   ! inverse M_3
   xyz2 = xyz3(1:2) * sqrt( 2.0*rs/(rs+abs(xyz3(3))) )
   
   ! inverse M_2
   qxy = sum(xyz2**2)
   
   special: if (dEq0(qxy)) then 
     Tinv = 0.0
   else special
     q2 = qxy + maxval(abs(xyz2))**2
     sq2 = sqrt(q2)
     q = (beta/r2/R1) * sqrt(q2*qxy/(q2-maxval(abs(xyz2))*sq2))
     tt = (minval(abs(xyz2))**2+maxval(abs(xyz2))*sq2)/r2/qxy 
     Tinv = q * sign(1.0,xyz2) * merge([ 1.0_pReal, acos(math_clip(tt,-1.0_pReal,1.0_pReal))/pi12], &
                                       [ acos(math_clip(tt,-1.0_pReal,1.0_pReal))/pi12, 1.0_pReal], &
                                       abs(xyz2(2)) <= abs(xyz2(1)))
   endif special
   
   ! inverse M_1
   xyz1 = [ Tinv(1), Tinv(2), sign(1.0,xyz3(3)) * rs / pref ] /sc
   
   ! reverst the coordinates back to the regular order according to the original pyramid number
   cube = xyz1(p)

 endif center

end function LambertBallToCube

!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief determine to which pyramid a point in a cubic grid belongs
!--------------------------------------------------------------------------
pure function GetPyramidOrder(xyz)

 implicit none
 real(pReal),intent(in),dimension(3) :: xyz
 integer(pInt),        dimension(3) :: GetPyramidOrder

 if      (((abs(xyz(1)) <=  xyz(3)).and.(abs(xyz(2)) <=  xyz(3))) .or. &
          ((abs(xyz(1)) <= -xyz(3)).and.(abs(xyz(2)) <= -xyz(3)))) then
   GetPyramidOrder = [1,2,3]
 else if (((abs(xyz(3)) <=  xyz(1)).and.(abs(xyz(2)) <=  xyz(1))) .or. &
          ((abs(xyz(3)) <= -xyz(1)).and.(abs(xyz(2)) <= -xyz(1)))) then
   GetPyramidOrder = [2,3,1]
 else if (((abs(xyz(1)) <=  xyz(2)).and.(abs(xyz(3)) <=  xyz(2))) .or. &
          ((abs(xyz(1)) <= -xyz(2)).and.(abs(xyz(3)) <= -xyz(2)))) then
   GetPyramidOrder = [3,1,2]
 else
   GetPyramidOrder = -1                                                                                  ! should be impossible, but might simplify debugging
 end if

end function GetPyramidOrder

end module Lambert
