! ###################################################################
! Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
! Modified      2017-2019, Martin Diehl/Max-Planck-Institut für Eisenforschung GmbH
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
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Mapping homochoric <-> cubochoric 
!
!> @details
!> D. Rosca, A. Morawiec, and M. De Graef. “A new method of constructing a grid 
!> in the space of 3D rotations and its applications to texture analysis”. 
!> Modeling and Simulations in Materials Science and Engineering 22, 075013 (2014).
!--------------------------------------------------------------------------
module Lambert
  use prec, only: &
    pReal
  use math, only: &
    PI

  implicit none
  private
  real(pReal), parameter, private :: &
    SPI  = sqrt(PI), &
    PREF = sqrt(6.0_pReal/PI), &
    A    = PI**(5.0_pReal/6.0_pReal)/6.0_pReal**(1.0_pReal/6.0_pReal), &
    AP   = PI**(2.0_pReal/3.0_pReal), &
    SC   = A/AP, &
    BETA = A/2.0_pReal, &
    R1   = (3.0_pReal*PI/4.0_pReal)**(1.0_pReal/3.0_pReal), &
    R2   = sqrt(2.0_pReal), &
    PI12 = PI/12.0_pReal, &
    PREK = R1 * 2.0_pReal**(1.0_pReal/4.0_pReal)/BETA
 
  public :: &
    LambertCubeToBall, &
    LambertBallToCube
  private :: &
   GetPyramidOrder

contains


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief map from 3D cubic grid to 3D ball 
!--------------------------------------------------------------------------
function LambertCubeToBall(cube) result(ball)
  use, intrinsic :: IEEE_ARITHMETIC
  use prec, only: &
    dEq0
 
  implicit none
  real(pReal), intent(in), dimension(3) :: cube
  real(pReal),             dimension(3) :: ball, LamXYZ, XYZ
  real(pReal),             dimension(2) :: T
  real(pReal)                           :: c, s, q
  real(pReal), parameter                :: eps = 1.0e-8_pReal
  integer,                 dimension(3) :: p
  integer,                 dimension(2) :: order
  
  if (maxval(abs(cube)) > AP/2.0+eps) then
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
      order = merge( [2,1], [1,2], abs(XYZ(2)) <= abs(XYZ(1)))                                      ! order of absolute values of XYZ
      q = PI12 *  XYZ(order(1))/XYZ(order(2))                                                       ! smaller by larger
      c = cos(q)
      s = sin(q)
      q = prek * XYZ(order(2))/ sqrt(R2-c)
      T = [ (R2*c - 1.0), R2 * s] * q
  
      ! transform to sphere grid (inverse Lambert)
      ! [note that there is no need to worry about dividing by zero, since XYZ(3) can not become zero]
      c = sum(T**2)
      s = Pi * c/(24.0*XYZ(3)**2)
      c = sPi * c / sqrt(24.0_pReal) / XYZ(3)
      q = sqrt( 1.0 - s )
      LamXYZ = [ T(order(2)) * q, T(order(1)) * q, pref * XYZ(3) - c ]
    endif special
  
    ! reverse the coordinates back to the regular order according to the original pyramid number
    ball = LamXYZ(p)
 
  endif center

end function LambertCubeToBall


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief map from 3D ball to 3D cubic grid  
!--------------------------------------------------------------------------
pure function LambertBallToCube(xyz) result(cube)
  use, intrinsic :: IEEE_ARITHMETIC, only:&
    IEEE_positive_inf, &
    IEEE_value
  use prec, only: &
    dEq0
   use math, only: &
    math_clip
 
  implicit none
  real(pReal), intent(in), dimension(3) :: xyz
  real(pReal),             dimension(3) :: cube, xyz1, xyz3
  real(pReal),             dimension(2) :: Tinv, xyz2
  real(pReal)                           :: rs, qxy, q2, sq2, q, tt
  integer,                 dimension(3) :: p
  
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
      Tinv = 0.0_pReal
    else special
      q2 = qxy + maxval(abs(xyz2))**2
      sq2 = sqrt(q2)
      q = (beta/R2/R1) * sqrt(q2*qxy/(q2-maxval(abs(xyz2))*sq2))
      tt = (minval(abs(xyz2))**2+maxval(abs(xyz2))*sq2)/R2/qxy 
      Tinv = q * sign(1.0_pReal,xyz2) * merge([ 1.0_pReal, acos(math_clip(tt,-1.0_pReal,1.0_pReal))/PI12], &
                                              [ acos(math_clip(tt,-1.0_pReal,1.0_pReal))/PI12, 1.0_pReal], &
                                               abs(xyz2(2)) <= abs(xyz2(1)))
    endif special
    
    ! inverse M_1
    xyz1 = [ Tinv(1), Tinv(2), sign(1.0_pReal,xyz3(3)) * rs / pref ] /sc
    
    ! reverst the coordinates back to the regular order according to the original pyramid number
    cube = xyz1(p)
 
  endif center

end function LambertBallToCube


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief determine to which pyramid a point in a cubic grid belongs
!--------------------------------------------------------------------------
pure function GetPyramidOrder(xyz)

  implicit none
  real(pReal),intent(in),dimension(3) :: xyz
  integer,               dimension(3) :: GetPyramidOrder
 
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
    GetPyramidOrder = -1                                                                            ! should be impossible, but might simplify debugging
  end if

end function GetPyramidOrder

end module Lambert
