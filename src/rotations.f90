! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
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

!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief rotation storage and conversion
!> @details: rotation is internally stored as quaternion. It cabe inialized from different 
!> represantations and also returns itself in different representations.
!
!  All methods and naming conventions based on Rowenhorst_etal2015
!    Convention 1: coordinate frames are right-handed
!    Convention 2: a rotation angle ω is taken to be positive for a counterclockwise rotation
!                  when viewing from the end point of the rotation axis towards the origin
!    Convention 3: rotations will be interpreted in the passive sense
!    Convention 4: Euler angle triplets are implemented using the Bunge convention,
!                  with the angular ranges as [0, 2π],[0, π],[0, 2π]
!    Convention 5: the rotation angle ω is limited to the interval [0, π]
!    Convention 6: P = -1
!---------------------------------------------------------------------------------------------------

module rotations
 use prec, only: &
   pReal
 use quaternions

 implicit none
 private
 type, public :: rotation
   type(quaternion), private :: q
   contains
     procedure, public :: asQuaternion
     procedure, public :: asEulerAngles
     procedure, public :: asAxisAnglePair
     procedure, public :: asRodriguesFrankVector
     procedure, public :: asRotationMatrix
     !------------------------------------------
     procedure, public :: fromRotationMatrix
     !------------------------------------------
     procedure, public :: rotVector
     procedure, public :: rotTensor
     procedure, public :: misorientation
 end type rotation


contains


!---------------------------------------------------------------------------------------------------
! Return rotation in different represenations
!---------------------------------------------------------------------------------------------------
function asQuaternion(self)

 implicit none
 class(rotation), intent(in) :: self
 real(pReal), dimension(4)   :: asQuaternion

 asQuaternion = [self%q%w, self%q%x, self%q%y, self%q%z]

end function asQuaternion
!---------------------------------------------------------------------------------------------------
function asEulerAngles(self)
  
 implicit none
 class(rotation), intent(in) :: self
 real(pReal), dimension(3)   :: asEulerAngles
  
 asEulerAngles = qu2eu(self%q)

end function asEulerAngles
!---------------------------------------------------------------------------------------------------
function asAxisAnglePair(self)

 implicit none
 class(rotation), intent(in) :: self
 real(pReal), dimension(4)   :: asAxisAnglePair

 asAxisAnglePair = qu2ax(self%q)

end function asAxisAnglePair
!---------------------------------------------------------------------------------------------------
function asRotationMatrix(self)
 
 implicit none
 class(rotation), intent(in) :: self
 real(pReal), dimension(3,3) :: asRotationMatrix

 asRotationMatrix = qu2om(self%q)

end function asRotationMatrix
!---------------------------------------------------------------------------------------------------
function asRodriguesFrankVector(self)

 implicit none
 class(rotation), intent(in) :: self
 real(pReal), dimension(4)   :: asRodriguesFrankVector

 asRodriguesFrankVector = qu2ro(self%q)
 
end function asRodriguesFrankVector
!---------------------------------------------------------------------------------------------------
function asHomochoric(self)

 implicit none
 class(rotation), intent(in) :: self
 real(pReal), dimension(3)   :: asHomochoric

 asHomochoric = qu2ho(self%q)

end function asHomochoric
 
!---------------------------------------------------------------------------------------------------
! Initialize rotation from different representations
!---------------------------------------------------------------------------------------------------
subroutine fromRotationMatrix(self,om)

 implicit none
 class(rotation), intent(out)            :: self
 real(pReal), dimension(3,3), intent(in) :: om

 self%q = om2qu(om)

end subroutine


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief rotate a vector passively (default) or actively
!> @details: rotation is based on unit quaternion or rotation matrix (fallback)
!---------------------------------------------------------------------------------------------------
function rotVector(self,v,active)
 use prec, only: &
   dEq
   
 implicit none
 real(pReal),                 dimension(3) :: rotVector
 class(rotation), intent(in)               :: self
 real(pReal),     intent(in), dimension(3) :: v
 logical,         intent(in), optional     :: active
 
 type(quaternion)                  :: q
  
 if (dEq(norm2(v),1.0_pReal,tol=1.0e-15_pReal)) then
   passive: if (merge(.not. active, .true., present(active))) then                                  ! ToDo: not save (PGI)
     q = self%q * (quaternion([0.0_pReal, v(1), v(2), v(3)]) * conjg(self%q) )
   else passive
     q = conjg(self%q) * (quaternion([0.0_pReal, v(1), v(2), v(3)]) * self%q )
   endif passive
   rotVector = [q%x,q%y,q%z]
 else
   passive2: if (merge(.not. active, .true., present(active))) then                                 ! ToDo: not save (PGI)
     rotVector = matmul(self%asRotationMatrix(),v)
   else passive2
     rotVector = matmul(transpose(self%asRotationMatrix()),v)
   endif passive2
 endif

end function rotVector


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief rotate a second rank tensor passively (default) or actively
!> @details: rotation is based on rotation matrix
!---------------------------------------------------------------------------------------------------
function rotTensor(self,m,active)
  
 implicit none
 real(pReal),                 dimension(3,3) :: rotTensor
 class(rotation), intent(in)                 :: self
 real(pReal),     intent(in), dimension(3,3) :: m
 logical,         intent(in), optional       :: active
  

 passive: if (merge(.not. active, .true., present(active))) then
   rotTensor = matmul(matmul(self%asRotationMatrix(),m),transpose(self%asRotationMatrix()))
 else passive
   rotTensor = matmul(matmul(transpose(self%asRotationMatrix()),m),self%asRotationMatrix())
 endif passive

end function rotTensor


!---------------------------------------------------------------------------------------------------
!> @brief misorientation
!---------------------------------------------------------------------------------------------------
function misorientation(self,other)
  
 implicit none
 type(rotation)              :: misorientation
 class(rotation), intent(in) :: self, other
 
 misorientation%q = conjg(self%q) * other%q                                                         !ToDo: this is the convention used in math

end function misorientation


!---------------------------------------------------------------------------------------------------
! The following routines convert between different representations
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Euler angles to orientation matrix
!---------------------------------------------------------------------------------------------------
pure function eu2om(eu) result(om)
 use prec, only: &
   dEq0
 
 implicit none
 real(pReal), intent(in), dimension(3)   :: eu
 real(pReal),             dimension(3,3) :: om
 
 real(pReal),             dimension(3)   :: c, s      

 c = cos(eu)
 s = sin(eu)

 om(1,1) =  c(1)*c(3)-s(1)*s(3)*c(2)
 om(1,2) =  s(1)*c(3)+c(1)*s(3)*c(2)
 om(1,3) =  s(3)*s(2)
 om(2,1) = -c(1)*s(3)-s(1)*c(3)*c(2)
 om(2,2) = -s(1)*s(3)+c(1)*c(3)*c(2)
 om(2,3) =  c(3)*s(2)
 om(3,1) =  s(1)*s(2)
 om(3,2) = -c(1)*s(2)
 om(3,3) =  c(2)

 where(dEq0(om)) om = 0.0_pReal

end function eu2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert euler to axis angle
!---------------------------------------------------------------------------------------------------
pure function eu2ax(eu) result(ax)
 use prec, only: &
   dEq0, &
   dEq
 use math, only: &
   PI

 implicit none
 real(pReal), intent(in), dimension(3) :: eu
 real(pReal),             dimension(4) :: ax
 
 real(pReal)                           :: t, delta, tau, alpha, sigma
 
 t     = tan(eu(2)*0.5)
 sigma = 0.5*(eu(1)+eu(3))
 delta = 0.5*(eu(1)-eu(3))
 tau   = sqrt(t**2+sin(sigma)**2)
 
 alpha = merge(PI, 2.0*atan(tau/cos(sigma)), dEq(sigma,PI*0.5_pReal,tol=1.0e-15_pReal))
 
 if (dEq0(alpha)) then                                                                              ! return a default identity axis-angle pair
   ax = [ 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal ]
 else
   ax(1:3) = -P/tau * [ t*cos(delta), t*sin(delta), sin(sigma) ]                                    ! passive axis-angle pair so a minus sign in front
   ax(4) = alpha
   if (alpha < 0.0) ax = -ax                                                                        ! ensure alpha is positive
 end if
 
end function eu2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Euler angles to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function eu2ro(eu) result(ro)
 use prec, only: &
   dEq0
 use, intrinsic :: IEEE_ARITHMETIC, only: &
   IEEE_value, &
   IEEE_positive_inf
 use math, only: &
   PI

 implicit none
 real(pReal), intent(in), dimension(3) :: eu
 real(pReal),             dimension(4) :: ro
 
 ro = eu2ax(eu)
 if (ro(4) >= PI) then
   ro(4) = IEEE_value(ro(4),IEEE_positive_inf)
 elseif(dEq0(ro(4))) then
   ro = [ 0.0_pReal, 0.0_pReal, P, 0.0_pReal ]
 else
   ro(4) = tan(ro(4)*0.5)
 end if
 
end function eu2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Euler angles to unit quaternion
!---------------------------------------------------------------------------------------------------
pure function eu2qu(eu) result(qu)

 implicit none
 real(pReal),      intent(in), dimension(3) :: eu
 type(quaternion)                           :: qu
 real(pReal),                  dimension(3) :: ee
 real(pReal)                                :: cPhi, sPhi

 ee = 0.5_pReal*eu
 
 cPhi = cos(ee(2))
 sPhi = sin(ee(2))

 qu = quaternion([   cPhi*cos(ee(1)+ee(3)), &
                  -P*sPhi*cos(ee(1)-ee(3)), &
                  -P*sPhi*sin(ee(1)-ee(3)), &
                  -P*cPhi*sin(ee(1)+ee(3))])
 if(qu%w < 0.0_pReal) qu = qu%homomorphed() 

end function eu2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief orientation matrix to Euler angles
!---------------------------------------------------------------------------------------------------
pure function om2eu(om) result(eu)
 use prec, only: &
   dEq
 use math, only: &
   PI

 implicit none
 real(pReal), intent(in), dimension(3,3) :: om
 real(pReal),             dimension(3)   :: eu
 real(pReal)                             :: zeta
 
 if (dEq(abs(om(3,3)),1.0_pReal,1.0e-15_pReal)) then
   eu = [ atan2( om(1,2),om(1,1)), 0.5*PI*(1-om(3,3)),0.0_pReal ]
 else 
   zeta = 1.0_pReal/sqrt(1.0_pReal-om(3,3)**2.0_pReal)
   eu = [atan2(om(3,1)*zeta,-om(3,2)*zeta), &
         acos(om(3,3)), &
         atan2(om(1,3)*zeta, om(2,3)*zeta)]
 end if
 where(eu<0.0_pReal) eu = mod(eu+2.0_pReal*PI,[2.0_pReal*PI,PI,2.0_pReal*PI])
 
end function om2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to orientation matrix
!---------------------------------------------------------------------------------------------------
pure function ax2om(ax) result(om)
 use prec, only: &
   pInt

 implicit none
 real(pReal), intent(in), dimension(4)   :: ax
 real(pReal),             dimension(3,3) :: om
 
 real(pReal)                             :: q, c, s, omc
 integer(pInt)                           :: i

 c = cos(ax(4))
 s = sin(ax(4))
 omc = 1.0-c

 forall(i=1:3) om(i,i) = ax(i)**2*omc + c

 q = omc*ax(1)*ax(2)
 om(1,2) = q + s*ax(3)
 om(2,1) = q - s*ax(3)
 
 q = omc*ax(2)*ax(3)
 om(2,3) = q + s*ax(1)
 om(3,2) = q - s*ax(1)
 
 q = omc*ax(3)*ax(1)
 om(3,1) = q + s*ax(2)
 om(1,3) = q - s*ax(2)

 if (P > 0.0) om = transpose(om)

end function ax2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to Euler angles
!---------------------------------------------------------------------------------------------------
pure function qu2eu(qu) result(eu)
 use prec, only: &
   dEq0
 use math, only: &
   PI

 implicit none
 type(quaternion), intent(in)              :: qu
 real(pReal),                 dimension(3) :: eu
 
 real(pReal)                               :: q12, q03, chi, chiInv
 
 q03 = qu%w**2+qu%z**2
 q12 = qu%x**2+qu%y**2
 chi = sqrt(q03*q12)
 
 degenerated: if (dEq0(chi)) then
   eu = merge([atan2(-P*2.0*qu%w*qu%z,qu%w**2-qu%z**2), 0.0_pReal, 0.0_pReal], &
              [atan2(2.0*qu%x*qu%y,qu%x**2-qu%y**2),    PI,        0.0_pReal], &
              dEq0(q12))
 else degenerated
   chiInv = 1.0/chi
   eu = [atan2((-P*qu%w*qu%y+qu%x*qu%z)*chi, (-P*qu%w*qu%x-qu%y*qu%z)*chi ), &
         atan2( 2.0*chi, q03-q12 ), &
         atan2(( P*qu%w*qu%y+qu%x*qu%z)*chi, (-P*qu%w*qu%x+qu%y*qu%z)*chi )]
 endif degenerated
 where(eu<0.0_pReal) eu = mod(eu+2.0_pReal*PI,[2.0_pReal*PI,PI,2.0_pReal*PI])

end function qu2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to homochoric
!---------------------------------------------------------------------------------------------------
pure function ax2ho(ax) result(ho)

 implicit none 
 real(pReal), intent(in), dimension(4) :: ax
 real(pReal),             dimension(3) :: ho
 
 real(pReal)                           :: f
 
 f = 0.75 * ( ax(4) - sin(ax(4)) )
 f = f**(1.0/3.0)
 ho = ax(1:3) * f

end function ax2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to axis angle pair
!---------------------------------------------------------------------------------------------------
pure function ho2ax(ho) result(ax)
 use prec, only: &
   pInt, &
   dEq0
 implicit none 
 real(pReal), intent(in), dimension(3) :: ho
 real(pReal),             dimension(4) :: ax
 
 integer(pInt)                         :: i
 real(pReal)                           :: hmag_squared, s, hm
 real(pReal), parameter, dimension(16) :: &
   tfit = [ 1.0000000000018852_pReal,      -0.5000000002194847_pReal, & 
           -0.024999992127593126_pReal,    -0.003928701544781374_pReal, & 
           -0.0008152701535450438_pReal,   -0.0002009500426119712_pReal, & 
           -0.00002397986776071756_pReal,  -0.00008202868926605841_pReal, & 
           +0.00012448715042090092_pReal,  -0.0001749114214822577_pReal, & 
           +0.0001703481934140054_pReal,   -0.00012062065004116828_pReal, & 
           +0.000059719705868660826_pReal, -0.00001980756723965647_pReal, & 
           +0.000003953714684212874_pReal, -0.00000036555001439719544_pReal ]
 
 ! normalize h and store the magnitude
 hmag_squared = sum(ho**2.0_pReal)
 if (dEq0(hmag_squared)) then
   ax = [ 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal ]
 else
   hm = hmag_squared
 
 ! convert the magnitude to the rotation angle
   s = tfit(1) + tfit(2) * hmag_squared
   do i=3,16
     hm = hm*hmag_squared
     s  = s + tfit(i) * hm
   end do
   ax = [ho/sqrt(hmag_squared), 2.0_pReal*acos(s)]
 end if

end function ho2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert orientation matrix to axis angle pair
!---------------------------------------------------------------------------------------------------
function om2ax(om) result(ax)
 use prec, only: &
   pInt, &
   dEq0, &
   cEq, &
   dNeq0
 use IO, only: &
   IO_error
 use math, only: &
   math_clip, &
   math_trace33

 implicit none
 real(pReal), intent(in)     :: om(3,3)
 real(pReal)                 :: ax(4)
 
 real(pReal)                 :: t
 real(pReal), dimension(3)   :: Wr, Wi
 real(pReal), dimension(10)  :: WORK
 real(pReal), dimension(3,3) :: VR, devNull, o
 integer(pInt)               :: INFO, LWORK, i
 
 external :: dgeev,sgeev
 
 o = om
 
 ! first get the rotation angle
 t = 0.5_pReal * (math_trace33(om) - 1.0)
 ax(4) = acos(math_clip(t,-1.0_pReal,1.0_pReal))
 
 if (dEq0(ax(4))) then
   ax(1:3) = [ 0.0, 0.0, 1.0 ]
 else
   ! set some initial LAPACK variables 
   INFO = 0
   ! first initialize the parameters for the LAPACK DGEEV routines
   LWORK = 20   
 
 ! call the eigenvalue solver
#if (FLOAT==8)
   call dgeev('N','V',3,o,3,Wr,Wi,devNull,3,VR,3,WORK,LWORK,INFO)
#elif (FLOAT==4)
   call sgeev('N','V',3,o,3,Wr,Wi,devNull,3,VR,3,WORK,LWORK,INFO)
#else
   NO SUITABLE PRECISION FOR REAL SELECTED, STOPPING COMPILATION
#endif
   if (INFO /= 0) call IO_error(0_pInt,ext_msg='Error in om2ax/(s/d)geev: (S/D)GEEV return not zero')
   i = maxloc(merge(1.0_pReal,0.0_pReal,cEq(cmplx(Wr,Wi,pReal),cmplx(1.0_pReal,0.0_pReal,pReal),tol=1.0e-14_pReal)),dim=1) ! poor substitute for findloc
   ax(1:3) = VR(1:3,i)
   where (                dNeq0([om(2,3)-om(3,2), om(3,1)-om(1,3), om(1,2)-om(2,1)])) &
     ax(1:3) = sign(ax(1:3),-P *[om(2,3)-om(3,2), om(3,1)-om(1,3), om(1,2)-om(2,1)])
 endif

end function om2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to axis angle pair
!---------------------------------------------------------------------------------------------------
pure function ro2ax(ro) result(ax)
 use, intrinsic :: IEEE_ARITHMETIC, only: &
   IEEE_is_finite
 use prec, only: &
   dEq0
 use math, only: &
   PI

 implicit none 
 real(pReal), intent(in), dimension(4) :: ro
 real(pReal),             dimension(4) :: ax
 
 real(pReal)                           :: ta, angle
 
 ta = ro(4)
 
 if (dEq0(ta))  then 
   ax = [ 0.0, 0.0, 1.0, 0.0 ]
 elseif (.not. IEEE_is_finite(ta)) then
   ax = [ ro(1), ro(2), ro(3), PI ]
 else
   angle = 2.0*atan(ta)
   ta = 1.0/norm2(ro(1:3))
   ax = [ ro(1)/ta, ro(2)/ta, ro(3)/ta, angle ]
 end if

end function ro2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function ax2ro(ax) result(ro)
 use, intrinsic :: IEEE_ARITHMETIC, only: &
   IEEE_value, &
   IEEE_positive_inf
 use prec, only: &
   dEq0
 use math, only: &
   PI

 implicit none 
 real(pReal), intent(in), dimension(4) :: ax
 real(pReal),             dimension(4) :: ro
 
 real(pReal), parameter                :: thr = 1.0E-7
 
 if (dEq0(ax(4))) then
   ro = [ 0.0_pReal, 0.0_pReal, P, 0.0_pReal ]
 else 
   ro(1:3) =  ax(1:3)
   ! we need to deal with the 180 degree case
   ro(4) = merge(IEEE_value(ro(4),IEEE_positive_inf),tan(ax(4)*0.5 ),abs(ax(4)-PI) < thr)
 end if

end function ax2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to quaternion
!---------------------------------------------------------------------------------------------------
pure function ax2qu(ax) result(qu)
 use prec, only: &
   dEq0
   
 implicit none
 real(pReal),      intent(in), dimension(4) :: ax
 type(quaternion)                           :: qu
 
 real(pReal)                                :: c, s


 if (dEq0(ax(4))) then
   qu = quaternion([ 1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal ])
 else
   c = cos(ax(4)*0.5)
   s = sin(ax(4)*0.5)
   qu = quaternion([ c, ax(1)*s, ax(2)*s, ax(3)*s ])
 end if

end function ax2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to homochoric
!---------------------------------------------------------------------------------------------------
pure function ro2ho(ro) result(ho)
 use, intrinsic :: IEEE_ARITHMETIC, only: &
   IEEE_is_finite
 use prec, only: &
   dEq0
 use math, only: &
   PI

 implicit none
 real(pReal),        intent(in), dimension(4) :: ro
 real(pReal),                    dimension(3) :: ho
 
 real(pReal)                                  :: f
 
 if (dEq0(norm2(ro(1:3)))) then
   ho = [ 0.0, 0.0, 0.0 ]
 else
   f = merge(2.0*atan(ro(4)) - sin(2.0*atan(ro(4))),PI, IEEE_is_finite(ro(4)))
   ho = ro(1:3) * (0.75_pReal*f)**(1.0/3.0)
 end if

end function ro2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to rotation matrix
!---------------------------------------------------------------------------------------------------
pure function qu2om(qu) result(om)

 implicit none
 type(quaternion), intent(in)                :: qu
 real(pReal),                 dimension(3,3) :: om
 
 real(pReal)                                 :: qq

 qq = qu%w**2-(qu%x**2 + qu%y**2 + qu%z**2)


 om(1,1) = qq+2.0*qu%x*qu%x
 om(2,2) = qq+2.0*qu%y*qu%y
 om(3,3) = qq+2.0*qu%z*qu%z

 om(1,2) = 2.0*(qu%x*qu%y-qu%w*qu%z)
 om(2,3) = 2.0*(qu%y*qu%z-qu%w*qu%x)
 om(3,1) = 2.0*(qu%z*qu%x-qu%w*qu%y)
 om(2,1) = 2.0*(qu%y*qu%x+qu%w*qu%z)
 om(3,2) = 2.0*(qu%z*qu%y+qu%w*qu%x)
 om(1,3) = 2.0*(qu%x*qu%z+qu%w*qu%y)

 if (P < 0.0) om = transpose(om)

end function qu2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to a unit quaternion
!---------------------------------------------------------------------------------------------------
function om2qu(om) result(qu)
 use prec, only: &
   dEq

 implicit none
 real(pReal),      intent(in), dimension(3,3) :: om
 type(quaternion)                             :: qu
 
 real(pReal),                  dimension(4)   :: qu_A
 real(pReal),                  dimension(4)   :: s

 s = [+om(1,1) +om(2,2) +om(3,3) +1.0_pReal, &
      +om(1,1) -om(2,2) -om(3,3) +1.0_pReal, &
      -om(1,1) +om(2,2) -om(3,3) +1.0_pReal, &
      -om(1,1) -om(2,2) +om(3,3) +1.0_pReal]

 qu_A = sqrt(max(s,0.0_pReal))*0.5_pReal*[1.0_pReal,P,P,P]
 qu_A = qu_A/norm2(qu_A)

 if(any(dEq(abs(qu_A),1.0_pReal,1.0e-15_pReal))) &
   where (.not.(dEq(abs(qu_A),1.0_pReal,1.0e-15_pReal))) qu_A = 0.0_pReal

 if (om(3,2) < om(2,3)) qu_A(2) = -qu_A(2)
 if (om(1,3) < om(3,1)) qu_A(3) = -qu_A(3)
 if (om(2,1) < om(1,2)) qu_A(4) = -qu_A(4)
 
 qu = quaternion(qu_A)
 !qu_A = om2ax(om)
 !if(any(qu_A(1:3) * [qu%x,qu%y,qu%z] < 0.0)) print*, 'sign error' 

end function om2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to axis angle pair
!---------------------------------------------------------------------------------------------------
pure function qu2ax(qu) result(ax)
 use prec, only: &
   dEq0, &
   dNeq0
 use math, only: &
   PI, &
   math_clip

 implicit none
 type(quaternion), intent(in)              :: qu
 real(pReal),                 dimension(4) :: ax
 
 real(pReal)                               :: omega, s

 omega = 2.0 * acos(math_clip(qu%w,0.0_pReal,1.0_pReal))
 ! if the angle equals zero, then we return the rotation axis as [001]
 if (dEq0(omega)) then
   ax = [ 0.0, 0.0, 1.0, 0.0 ]
 elseif (dNeq0(qu%w)) then
   s =  sign(1.0_pReal,qu%w)/sqrt(qu%x**2+qu%y**2+qu%z**2)
   ax = [ qu%x*s, qu%y*s, qu%z*s, omega ]
 else
   ax = [ qu%x, qu%y, qu%z, PI ]
 end if

end function qu2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function qu2ro(qu) result(ro)
 use, intrinsic :: IEEE_ARITHMETIC, only: &
   IEEE_value, &
   IEEE_positive_inf
 use prec, only: &
   dEq0
 
 type(quaternion), intent(in)              :: qu
 real(pReal),                 dimension(4) :: ro
 
 real(pReal)                               :: s
 real(pReal), parameter                    :: thr = 1.0e-8_pReal
 
 if (qu%w < thr) then
   ro = [qu%x, qu%y, qu%z, IEEE_value(ro(4),IEEE_positive_inf)]
 else
   s = norm2([qu%x,qu%y,qu%z])
   ro = merge ( [ 0.0_pReal, 0.0_pReal, P, 0.0_pReal], &
                [ qu%x/s,  qu%y/s,  qu%z/s, tan(acos(qu%w))], &
                s < thr)                                                                            !ToDo: not save (PGI compiler)
 end if

end function qu2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to homochoric
!---------------------------------------------------------------------------------------------------
pure function qu2ho(qu) result(ho)
 use prec, only: &
   dEq0
 use math, only: &
   math_clip

 implicit none
 type(quaternion), intent(in)               :: qu
 real(pReal),                  dimension(3) :: ho
 
 real(pReal)                                :: omega, f

 omega = 2.0 * acos(math_clip(qu%w,0.0_pReal,1.0_pReal))
 
 if (dEq0(omega)) then
   ho = [ 0.0, 0.0, 0.0 ]
 else
   ho = [qu%x, qu%y, qu%z]
   f  = 0.75 * ( omega - sin(omega) )
   ho = ho/norm2(ho)* f**(1.0/3.0)
 end if

end function qu2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to cubochoric
!---------------------------------------------------------------------------------------------------
function ho2cu(ho) result(cu)
 use Lambert, only: &
   LambertBallToCube

 implicit none
 real(pReal), intent(in), dimension(3) :: ho
 real(pReal),             dimension(3) :: cu

 cu = LambertBallToCube(ho)

end function ho2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to homochoric
!---------------------------------------------------------------------------------------------------
function cu2ho(cu) result(ho)
 use Lambert, only: &
   LambertCubeToBall

 implicit none
 real(pReal), intent(in), dimension(3) :: cu
 real(pReal),             dimension(3) :: ho

 ho = LambertCubeToBall(cu)

end function cu2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to Euler angles
!---------------------------------------------------------------------------------------------------
pure function ro2eu(ro) result(eu)

 implicit none
 real(pReal), intent(in), dimension(4) :: ro
 real(pReal),             dimension(3) :: eu
        
 eu = om2eu(ro2om(ro))

end function ro2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Euler angles to homochoric
!---------------------------------------------------------------------------------------------------
pure function eu2ho(eu) result(ho)

 implicit none
 real(pReal), intent(in), dimension(3) :: eu
 real(pReal),             dimension(3) :: ho

 ho = ax2ho(eu2ax(eu))

end function eu2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function om2ro(om) result(ro)

 implicit none
 real(pReal), intent(in), dimension(3,3) :: om
 real(pReal),             dimension(4)   :: ro

 ro = eu2ro(om2eu(om))

end function om2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to homochoric
!---------------------------------------------------------------------------------------------------
function om2ho(om) result(ho)

 implicit none
 real(pReal), intent(in), dimension(3,3) :: om
 real(pReal),             dimension(3)   :: ho

 ho = ax2ho(om2ax(om))

end function om2ho 


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to Euler angles
!---------------------------------------------------------------------------------------------------
pure function ax2eu(ax) result(eu)

 implicit none
 real(pReal), intent(in), dimension(4) :: ax
 real(pReal),             dimension(3) :: eu

 eu = om2eu(ax2om(ax))

end function ax2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to rotation matrix
!---------------------------------------------------------------------------------------------------
pure function ro2om(ro) result(om)

 implicit none
 real(pReal), intent(in), dimension(4)   :: ro
 real(pReal),             dimension(3,3) :: om

 om = ax2om(ro2ax(ro))

end function ro2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to unit quaternion
!---------------------------------------------------------------------------------------------------
pure function ro2qu(ro) result(qu)

 implicit none
 real(pReal),      intent(in), dimension(4) :: ro
 type(quaternion)                           :: qu
 
 qu = ax2qu(ro2ax(ro))

end function ro2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to Euler angles
!---------------------------------------------------------------------------------------------------
pure function ho2eu(ho) result(eu)

 implicit none
 real(pReal), intent(in), dimension(3) :: ho
 real(pReal),             dimension(3) :: eu

 eu = ax2eu(ho2ax(ho))

end function ho2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to rotation matrix
!---------------------------------------------------------------------------------------------------
pure function ho2om(ho) result(om)

 implicit none
 real(pReal), intent(in), dimension(3)   :: ho
 real(pReal),             dimension(3,3) :: om

 om = ax2om(ho2ax(ho))

end function ho2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function ho2ro(ho) result(ro)

 implicit none
 real(pReal), intent(in), dimension(3) :: ho
 real(pReal),             dimension(4) :: ro


 ro = ax2ro(ho2ax(ho))

end function ho2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to unit quaternion
!---------------------------------------------------------------------------------------------------
pure function ho2qu(ho) result(qu)

 implicit none
 real(pReal),      intent(in), dimension(3) :: ho
 type(quaternion)                           :: qu

 qu = ax2qu(ho2ax(ho))

end function ho2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Euler angles to cubochoric
!---------------------------------------------------------------------------------------------------
function eu2cu(eu) result(cu)

 implicit none
 real(pReal), intent(in), dimension(3) :: eu
 real(pReal),             dimension(3) :: cu

 cu = ho2cu(eu2ho(eu))

end function eu2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to cubochoric
!---------------------------------------------------------------------------------------------------
function om2cu(om) result(cu)

 implicit none
 real(pReal), intent(in), dimension(3,3) :: om
 real(pReal),             dimension(3)   :: cu

 cu = ho2cu(om2ho(om))

end function om2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to cubochoric
!---------------------------------------------------------------------------------------------------
function ax2cu(ax) result(cu)

 implicit none
 real(pReal), intent(in), dimension(4) :: ax
 real(pReal),             dimension(3) :: cu

 cu  = ho2cu(ax2ho(ax))

end function ax2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to cubochoric
!---------------------------------------------------------------------------------------------------
function ro2cu(ro) result(cu)

 implicit none
 real(pReal), intent(in), dimension(4) :: ro
 real(pReal),             dimension(3) :: cu

 cu = ho2cu(ro2ho(ro))

end function ro2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to cubochoric
!---------------------------------------------------------------------------------------------------
function qu2cu(qu) result(cu)
 
 implicit none
 type(quaternion), intent(in)              :: qu
 real(pReal),                 dimension(3) :: cu

 cu = ho2cu(qu2ho(qu))

end function qu2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to Euler angles
!---------------------------------------------------------------------------------------------------
function cu2eu(cu) result(eu)

 implicit none
 real(pReal), intent(in), dimension(3) :: cu
 real(pReal),             dimension(3) :: eu

 eu = ho2eu(cu2ho(cu))

end function cu2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to rotation matrix
!---------------------------------------------------------------------------------------------------
function cu2om(cu) result(om)

 implicit none
 real(pReal), intent(in), dimension(3)   :: cu
 real(pReal),             dimension(3,3) :: om

 om = ho2om(cu2ho(cu))

end function cu2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to axis angle pair
!---------------------------------------------------------------------------------------------------
function cu2ax(cu) result(ax)

 implicit none
 real(pReal), intent(in), dimension(3) :: cu
 real(pReal),             dimension(4) :: ax

 ax = ho2ax(cu2ho(cu))

end function cu2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to Rodrigues vector
!---------------------------------------------------------------------------------------------------
function cu2ro(cu) result(ro)

 implicit none
 real(pReal), intent(in), dimension(3) :: cu
 real(pReal),             dimension(4) :: ro

 ro = ho2ro(cu2ho(cu))

end function cu2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to unit quaternion
!---------------------------------------------------------------------------------------------------
function cu2qu(cu) result(qu)

 implicit none
 real(pReal), intent(in), dimension(3) :: cu
 type(quaternion)                      :: qu

 qu = ho2qu(cu2ho(cu))

end function cu2qu

end module rotations
