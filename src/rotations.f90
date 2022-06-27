! ###################################################################
! Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
! Modified      2017-2020, Martin Diehl/Max-Planck-Institut für Eisenforschung GmbH
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

!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief rotation storage and conversion
!> @details: rotation is internally stored as quaternion. It can be inialized from different
!> representations and also returns itself in different representations.
!
!  All methods and naming conventions based on Rowenhorst et al. 2015
!    Convention 1: coordinate frames are right-handed
!    Convention 2: a rotation angle ω is taken to be positive for a counterclockwise rotation
!                  when viewing from the end point of the rotation axis towards the origin
!    Convention 3: rotations will be interpreted in the passive sense
!    Convention 4: Euler angle triplets are implemented using the Bunge convention,
!                  with the angular ranges as [0, 2π],[0, π],[0, 2π]
!    Convention 5: the rotation angle ω is limited to the interval [0, π]
!    Convention 6: the real part of a quaternion is positive, Re(q) > 0
!    Convention 7: P = -1
!--------------------------------------------------------------------------------------------------

module rotations
  use IO
  use math

  implicit none
  private

  real(pReal), parameter :: P = -1.0_pReal                                                          !< parameter for orientation conversion.

  type, public :: tRotation
    real(pReal), dimension(4) :: q
    contains
      procedure, public :: asQuaternion
      procedure, public :: asEulers
      procedure, public :: asAxisAngle
      procedure, public :: asMatrix
      !------------------------------------------
      procedure, public :: fromQuaternion
      procedure, public :: fromEulers
      procedure, public :: fromAxisAngle
      procedure, public :: fromMatrix
      !------------------------------------------
      procedure, private :: rotRot__
      generic,   public  :: operator(*) => rotRot__
      generic,   public  :: rotate => rotVector,rotTensor2,rotTensor4
      procedure, public  :: rotVector
      procedure, public  :: rotTensor2
      procedure, public  :: rotTensor4
      procedure, public  :: rotStiffness
      procedure, public  :: misorientation
      procedure, public  :: standardize
  end type tRotation

  real(pReal), parameter :: &
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
    rotations_init, &
    eu2om

contains

!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine rotations_init

  print'(/,1x,a)', '<<<+-  rotations init  -+>>>'; flush(IO_STDOUT)

  print'(/,1x,a)', 'D. Rowenhorst et al., Modelling and Simulation in Materials Science and Engineering 23:083501, 2015'
  print'(  1x,a)', 'https://doi.org/10.1088/0965-0393/23/8/083501'

  call selfTest

end subroutine rotations_init


!--------------------------------------------------------------------------------------------------
! Return rotation in different representations.
!--------------------------------------------------------------------------------------------------
pure function asQuaternion(self)

  class(tRotation), intent(in) :: self
  real(pReal), dimension(4)    :: asQuaternion


  asQuaternion = self%q

end function asQuaternion
!--------------------------------------------------------------------------------------------------
pure function asEulers(self)

  class(tRotation), intent(in) :: self
  real(pReal), dimension(3)    :: asEulers


  asEulers = qu2eu(self%q)

end function asEulers
!--------------------------------------------------------------------------------------------------
pure function asAxisAngle(self)

  class(tRotation), intent(in) :: self
  real(pReal), dimension(4)    :: asAxisAngle


  asAxisAngle = qu2ax(self%q)

end function asAxisAngle
!--------------------------------------------------------------------------------------------------
pure function asMatrix(self)

  class(tRotation), intent(in) :: self
  real(pReal), dimension(3,3)  :: asMatrix


  asMatrix = qu2om(self%q)

end function asMatrix

!--------------------------------------------------------------------------------------------------
! Initialize rotation from different representations.
!--------------------------------------------------------------------------------------------------
subroutine fromQuaternion(self,qu)

  class(tRotation), intent(out)         :: self
  real(pReal), dimension(4), intent(in) :: qu


  if (dNeq(norm2(qu),1.0_pReal,1.0e-8_pReal)) call IO_error(402,ext_msg='fromQuaternion')

  self%q = qu

end subroutine fromQuaternion
!--------------------------------------------------------------------------------------------------
subroutine fromEulers(self,eu,degrees)

  class(tRotation), intent(out)         :: self
  real(pReal), dimension(3), intent(in) :: eu
  logical, intent(in), optional         :: degrees

  real(pReal), dimension(3)             :: Eulers


  if (.not. present(degrees)) then
    Eulers = eu
  else
    Eulers = merge(eu*INRAD,eu,degrees)
  end if

  if (any(Eulers<0.0_pReal) .or. any(Eulers>TAU) .or. Eulers(2) > PI) &
    call IO_error(402,ext_msg='fromEulers')

  self%q = eu2qu(Eulers)

end subroutine fromEulers
!--------------------------------------------------------------------------------------------------
subroutine fromAxisAngle(self,ax,degrees,P)

  class(tRotation), intent(out)         :: self
  real(pReal), dimension(4), intent(in) :: ax
  logical, intent(in), optional         :: degrees
  integer, intent(in), optional         :: P

  real(pReal)                           :: angle
  real(pReal),dimension(3)              :: axis


  if (.not. present(degrees)) then
    angle = ax(4)
  else
    angle = merge(ax(4)*INRAD,ax(4),degrees)
  end if

  if (.not. present(P)) then
    axis = ax(1:3)
  else
    axis = ax(1:3) * merge(-1.0_pReal,1.0_pReal,P == 1)
    if(abs(P) /= 1) call IO_error(402,ext_msg='fromAxisAngle (P)')
  end if

  if(dNeq(norm2(axis),1.0_pReal) .or. angle < 0.0_pReal .or. angle > PI) &
    call IO_error(402,ext_msg='fromAxisAngle')

  self%q = ax2qu([axis,angle])

end subroutine fromAxisAngle
!--------------------------------------------------------------------------------------------------
subroutine fromMatrix(self,om)

  class(tRotation), intent(out)           :: self
  real(pReal), dimension(3,3), intent(in) :: om


  if (dNeq(math_det33(om),1.0_pReal,tol=1.0e-5_pReal)) &
    call IO_error(402,ext_msg='fromMatrix')

  self%q = om2qu(om)

end subroutine fromMatrix
!--------------------------------------------------------------------------------------------------


!--------------------------------------------------------------------------------------------------
!> @brief: Compose rotations.
!--------------------------------------------------------------------------------------------------
pure elemental function rotRot__(self,R) result(rRot)

  type(tRotation)              :: rRot
  class(tRotation), intent(in) :: self,R


  rRot = tRotation(multiplyQuaternion(self%q,R%q))
  call rRot%standardize()

end function rotRot__


!--------------------------------------------------------------------------------------------------
!> @brief Convert to quaternion representation with positive q(1).
!--------------------------------------------------------------------------------------------------
pure elemental subroutine standardize(self)

  class(tRotation), intent(inout) :: self


  if (sign(1.0_pReal,self%q(1)) < 0.0_pReal) self%q = - self%q

end subroutine standardize


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Rotate a vector passively (default) or actively.
!--------------------------------------------------------------------------------------------------
pure function rotVector(self,v,active) result(vRot)

  real(pReal),                 dimension(3) :: vRot
  class(tRotation), intent(in)              :: self
  real(pReal),     intent(in), dimension(3) :: v
  logical,         intent(in), optional     :: active

  real(pReal), dimension(4) :: v_normed, q
  logical                   :: passive


  if (present(active)) then
    passive = .not. active
  else
    passive = .true.
  end if

  if (dEq0(norm2(v))) then
    vRot = v
  else
    v_normed = [0.0_pReal,v]/norm2(v)
    q = merge(multiplyQuaternion(self%q, multiplyQuaternion(v_normed, conjugateQuaternion(self%q))), &
              multiplyQuaternion(conjugateQuaternion(self%q), multiplyQuaternion(v_normed, self%q)), &
              passive)
    vRot = q(2:4)*norm2(v)
  end if

end function rotVector


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Rotate a rank-2 tensor passively (default) or actively.
!> @details: Rotation is based on rotation matrix
!--------------------------------------------------------------------------------------------------
pure function rotTensor2(self,T,active) result(tRot)

  real(pReal),                 dimension(3,3) :: tRot
  class(tRotation), intent(in)                :: self
  real(pReal),     intent(in), dimension(3,3) :: T
  logical,         intent(in), optional       :: active

  logical           :: passive


  if (present(active)) then
    passive = .not. active
  else
    passive = .true.
  end if

  tRot = merge(matmul(matmul(self%asMatrix(),T),transpose(self%asMatrix())), &
               matmul(matmul(transpose(self%asMatrix()),T),self%asMatrix()), &
               passive)

end function rotTensor2


!--------------------------------------------------------------------------------------------------
!> @brief Rotate a rank-4 tensor passively (default) or actively.
!> @details: rotation is based on rotation matrix
!! ToDo: Need to check active/passive !!!
!--------------------------------------------------------------------------------------------------
pure function rotTensor4(self,T,active) result(tRot)

  real(pReal),                 dimension(3,3,3,3) :: tRot
  class(tRotation), intent(in)                    :: self
  real(pReal),     intent(in), dimension(3,3,3,3) :: T
  logical,         intent(in), optional           :: active

  real(pReal), dimension(3,3) :: R
  integer :: i,j,k,l,m,n,o,p


  if (present(active)) then
    R = merge(transpose(self%asMatrix()),self%asMatrix(),active)
  else
    R = self%asMatrix()
  end if

  tRot = 0.0_pReal
  do i = 1,3;do j = 1,3;do k = 1,3;do l = 1,3
  do m = 1,3;do n = 1,3;do o = 1,3;do p = 1,3
    tRot(i,j,k,l) = tRot(i,j,k,l) &
                  + R(i,m) * R(j,n) * R(k,o) * R(l,p) * T(m,n,o,p)
  end do; end do; end do; end do; end do; end do; end do; end do

end function rotTensor4


!--------------------------------------------------------------------------------------------------
!> @brief Rotate a rank-4 stiffness tensor in Voigt 6x6 notation passively (default) or actively.
!> @details: https://scicomp.stackexchange.com/questions/35600
!! ToDo: Need to check active/passive !!!
!--------------------------------------------------------------------------------------------------
pure function rotStiffness(self,C,active) result(cRot)

  real(pReal),                 dimension(6,6) :: cRot
  class(tRotation), intent(in)                :: self
  real(pReal),     intent(in), dimension(6,6) :: C
  logical,         intent(in), optional       :: active

  real(pReal), dimension(3,3) :: R
  real(pReal), dimension(6,6) :: M


  if (present(active)) then
    R = merge(transpose(self%asMatrix()),self%asMatrix(),active)
  else
    R = self%asMatrix()
  end if

  M = reshape([R(1,1)**2,                   R(2,1)**2,                   R(3,1)**2, &
               R(2,1)*R(3,1),               R(1,1)*R(3,1),               R(1,1)*R(2,1), &
               R(1,2)**2,                   R(2,2)**2,                   R(3,2)**2, &
               R(2,2)*R(3,2),               R(1,2)*R(3,2),               R(1,2)*R(2,2), &
               R(1,3)**2,                   R(2,3)**2,                   R(3,3)**2, &
               R(2,3)*R(3,3),               R(1,3)*R(3,3),               R(1,3)*R(2,3), &
               2.0_pReal*R(1,2)*R(1,3),     2.0_pReal*R(2,2)*R(2,3),     2.0_pReal*R(3,2)*R(3,3), &
               R(2,2)*R(3,3)+R(2,3)*R(3,2), R(1,2)*R(3,3)+R(1,3)*R(3,2), R(1,2)*R(2,3)+R(1,3)*R(2,2), &
               2.0_pReal*R(1,3)*R(1,1),     2.0_pReal*R(2,3)*R(2,1),     2.0_pReal*R(3,3)*R(3,1), &
               R(2,3)*R(3,1)+R(2,1)*R(3,3), R(1,3)*R(3,1)+R(1,1)*R(3,3), R(1,3)*R(2,1)+R(1,1)*R(2,3), &
               2.0_pReal*R(1,1)*R(1,2),     2.0_pReal*R(2,1)*R(2,2),     2.0_pReal*R(3,1)*R(3,2), &
               R(2,1)*R(3,2)+R(2,2)*R(3,1), R(1,1)*R(3,2)+R(1,2)*R(3,1), R(1,1)*R(2,2)+R(1,2)*R(2,1)],[6,6])

  cRot = matmul(M,matmul(C,transpose(M)))

end function rotStiffness


!--------------------------------------------------------------------------------------------------
!> @brief Calculate misorientation.
!--------------------------------------------------------------------------------------------------
pure elemental function misorientation(self,other)

  type(tRotation)              :: misorientation
  class(tRotation), intent(in) :: self, other


  misorientation%q = multiplyQuaternion(other%q, conjugateQuaternion(self%q))

end function misorientation


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert unit quaternion to rotation matrix.
!--------------------------------------------------------------------------------------------------
pure function qu2om(qu) result(om)

  real(pReal), intent(in), dimension(4)   :: qu
  real(pReal),             dimension(3,3) :: om

  real(pReal)                             :: qq


  qq = qu(1)**2-sum(qu(2:4)**2)

  om(1,1) = qq+2.0_pReal*qu(2)**2
  om(2,2) = qq+2.0_pReal*qu(3)**2
  om(3,3) = qq+2.0_pReal*qu(4)**2

  om(1,2) = 2.0_pReal*(qu(2)*qu(3)-qu(1)*qu(4))
  om(2,3) = 2.0_pReal*(qu(3)*qu(4)-qu(1)*qu(2))
  om(3,1) = 2.0_pReal*(qu(4)*qu(2)-qu(1)*qu(3))
  om(2,1) = 2.0_pReal*(qu(3)*qu(2)+qu(1)*qu(4))
  om(3,2) = 2.0_pReal*(qu(4)*qu(3)+qu(1)*qu(2))
  om(1,3) = 2.0_pReal*(qu(2)*qu(4)+qu(1)*qu(3))

  if (sign(1.0_pReal,P) < 0.0_pReal) om = transpose(om)
  om = om/math_det33(om)**(1.0_pReal/3.0_pReal)

end function qu2om


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert unit quaternion to Bunge Euler angles.
!--------------------------------------------------------------------------------------------------
pure function qu2eu(qu) result(eu)

  real(pReal), intent(in), dimension(4) :: qu
  real(pReal),             dimension(3) :: eu

  real(pReal)                           :: q12, q03, chi


  q03 = qu(1)**2+qu(4)**2
  q12 = qu(2)**2+qu(3)**2
  chi = sqrt(q03*q12)

  degenerated: if (dEq0(q12)) then
    eu = [atan2(-P*2.0_pReal*qu(1)*qu(4),qu(1)**2-qu(4)**2), 0.0_pReal, 0.0_pReal]
  elseif          (dEq0(q03)) then
    eu = [atan2(   2.0_pReal*qu(2)*qu(3),qu(2)**2-qu(3)**2), PI,        0.0_pReal]
  else degenerated
    eu = [atan2((-P*qu(1)*qu(3)+qu(2)*qu(4))*chi, (-P*qu(1)*qu(2)-qu(3)*qu(4))*chi ), &
          atan2( 2.0_pReal*chi, q03-q12 ), &
          atan2(( P*qu(1)*qu(3)+qu(2)*qu(4))*chi, (-P*qu(1)*qu(2)+qu(3)*qu(4))*chi )]
  end if degenerated
  where(sign(1.0_pReal,eu)<0.0_pReal) eu = mod(eu+TAU,[TAU,PI,TAU])

end function qu2eu


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert unit quaternion to axis-angle pair.
!--------------------------------------------------------------------------------------------------
pure function qu2ax(qu) result(ax)

  real(pReal), intent(in), dimension(4) :: qu
  real(pReal),             dimension(4) :: ax

  real(pReal)                           :: omega, s


  if (dEq0(sum(qu(2:4)**2))) then
    ax = [ 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal ]                                             ! axis = [001]
  elseif (dNeq0(qu(1))) then
    s =  sign(1.0_pReal,qu(1))/norm2(qu(2:4))
    omega = 2.0_pReal * acos(math_clip(qu(1),-1.0_pReal,1.0_pReal))
    ax = [ qu(2)*s, qu(3)*s, qu(4)*s, omega ]
  else
    ax = [ qu(2),   qu(3),   qu(4),   PI ]
  end if

end function qu2ax


!--------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Convert rotation matrix to unit quaternion.
!> @details the original formulation (direct conversion) had (numerical?) issues
!--------------------------------------------------------------------------------------------------
pure function om2qu(om) result(qu)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(4)   :: qu

  real(pReal) :: trace,s
  trace = math_trace33(om)


  if(trace > 0.0_pReal) then
    s = 0.5_pReal / sqrt(trace+1.0_pReal)
    qu = [0.25_pReal/s, (om(3,2)-om(2,3))*s,(om(1,3)-om(3,1))*s,(om(2,1)-om(1,2))*s]
  else
      if( om(1,1) > om(2,2) .and. om(1,1) > om(3,3) ) then
          s = 2.0_pReal * sqrt( 1.0_pReal + om(1,1) - om(2,2) - om(3,3))
          qu = [ (om(3,2) - om(2,3)) /s,0.25_pReal * s,(om(1,2) + om(2,1)) / s,(om(1,3) + om(3,1)) / s]
      elseif (om(2,2) > om(3,3)) then
          s = 2.0_pReal * sqrt( 1.0_pReal + om(2,2) - om(1,1) - om(3,3))
          qu = [ (om(1,3) - om(3,1)) /s,(om(1,2) + om(2,1)) / s,0.25_pReal * s,(om(2,3) + om(3,2)) / s]
      else
          s = 2.0_pReal * sqrt( 1.0_pReal + om(3,3) - om(1,1) - om(2,2) )
          qu = [ (om(2,1) - om(1,2)) /s,(om(1,3) + om(3,1)) / s,(om(2,3) + om(3,2)) / s,0.25_pReal * s]
      end if
  end if
  if(sign(1.0_pReal,qu(1))<0.0_pReal) qu =-1.0_pReal * qu
  qu(2:4) = merge(qu(2:4),qu(2:4)*P,dEq0(qu(2:4)))
  qu = qu/norm2(qu)

end function om2qu


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert orientation matrix to Bunge Euler angles.
!> @details Two step check for special cases to avoid invalid operations (not needed for python)
!--------------------------------------------------------------------------------------------------
pure function om2eu(om) result(eu)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(3)   :: eu
  real(pReal)                             :: zeta


  if    (dNeq(abs(om(3,3)),1.0_pReal,1.e-8_pReal)) then
    zeta = 1.0_pReal/sqrt(math_clip(1.0_pReal-om(3,3)**2,1e-64_pReal,1.0_pReal))
    eu = [atan2(om(3,1)*zeta,-om(3,2)*zeta), &
          acos(math_clip(om(3,3),-1.0_pReal,1.0_pReal)), &
          atan2(om(1,3)*zeta, om(2,3)*zeta)]
  else
    eu = [atan2(om(1,2),om(1,1)), 0.5_pReal*PI*(1.0_pReal-om(3,3)),0.0_pReal ]
  end if
  where(abs(eu) < 1.e-8_pReal) eu = 0.0_pReal
  where(sign(1.0_pReal,eu)<0.0_pReal) eu = mod(eu+TAU,[TAU,PI,TAU])

end function om2eu


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert orientation matrix to axis-angle pair.
!--------------------------------------------------------------------------------------------------
function om2ax(om) result(ax)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(4)   :: ax

  real(pReal)                      :: t
  real(pReal), dimension(3)        :: Wr, Wi
  real(pReal), dimension((64+2)*3) :: work
  real(pReal), dimension(3,3)      :: VR, devNull, om_
  integer                          :: ierr, i


  om_ = om

  ! first get the rotation angle
  t = 0.5_pReal * (math_trace33(om) - 1.0_pReal)
  ax(4) = acos(math_clip(t,-1.0_pReal,1.0_pReal))

  if (dEq0(ax(4))) then
    ax(1:3) = [ 0.0_pReal, 0.0_pReal, 1.0_pReal ]
  else
    call dgeev('N','V',3,om_,3,Wr,Wi,devNull,3,VR,3,work,size(work,1),ierr)
    if (ierr /= 0) error stop 'LAPACK error'
    i = findloc(cEq(cmplx(Wr,Wi,pReal),cmplx(1.0_pReal,0.0_pReal,pReal),tol=1.0e-14_pReal),.true.,dim=1) !find eigenvalue (1,0)
    if (i == 0) error stop 'om2ax conversion failed'
    ax(1:3) = VR(1:3,i)
    where (                dNeq0([om(2,3)-om(3,2), om(3,1)-om(1,3), om(1,2)-om(2,1)])) &
      ax(1:3) = sign(ax(1:3),-P *[om(2,3)-om(3,2), om(3,1)-om(1,3), om(1,2)-om(2,1)])
  end if

end function om2ax


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert Bunge Euler angles to unit quaternion.
!--------------------------------------------------------------------------------------------------
pure function eu2qu(eu) result(qu)

  real(pReal), intent(in), dimension(3) :: eu
  real(pReal),             dimension(4) :: qu
  real(pReal),             dimension(3) :: ee
  real(pReal)                           :: cPhi, sPhi


  ee = 0.5_pReal*eu

  cPhi = cos(ee(2))
  sPhi = sin(ee(2))

  qu = [   cPhi*cos(ee(1)+ee(3)), &
        -P*sPhi*cos(ee(1)-ee(3)), &
        -P*sPhi*sin(ee(1)-ee(3)), &
        -P*cPhi*sin(ee(1)+ee(3))]
  if(sign(1.0_pReal,qu(1)) < 0.0_pReal) qu = qu * (-1.0_pReal)

end function eu2qu


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert Euler angles to orientation matrix.
!--------------------------------------------------------------------------------------------------
pure function eu2om(eu) result(om)

  real(pReal), intent(in), dimension(3)   :: eu
  real(pReal),             dimension(3,3) :: om

  real(pReal),             dimension(3)   :: c, s


  c = cos(eu)
  s = sin(eu)

  om(1,1) =  c(1)*c(3)-s(1)*s(3)*c(2)
  om(2,1) = -c(1)*s(3)-s(1)*c(3)*c(2)
  om(3,1) =  s(1)*s(2)
  om(1,2) =  s(1)*c(3)+c(1)*s(3)*c(2)
  om(2,2) = -s(1)*s(3)+c(1)*c(3)*c(2)
  om(3,2) = -c(1)*s(2)
  om(1,3) =  s(3)*s(2)
  om(2,3) =  c(3)*s(2)
  om(3,3) =  c(2)

  where(abs(om)<1.0e-12_pReal) om = 0.0_pReal

end function eu2om


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert Bunge Euler angles to axis-angle pair.
!--------------------------------------------------------------------------------------------------
pure function eu2ax(eu) result(ax)

  real(pReal), intent(in), dimension(3) :: eu
  real(pReal),             dimension(4) :: ax

  real(pReal)                           :: t, delta, tau, alpha, sigma


  t     = tan(eu(2)*0.5_pReal)
  sigma = 0.5_pReal*(eu(1)+eu(3))
  delta = 0.5_pReal*(eu(1)-eu(3))
  tau   = sqrt(t**2+sin(sigma)**2)

  alpha = merge(PI, 2.0_pReal*atan(tau/cos(sigma)), dEq(sigma,PI*0.5_pReal,tol=1.0e-15_pReal))

  if (dEq0(alpha)) then                                                                             ! return a default identity axis-angle pair
    ax = [ 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal ]
  else
    ax(1:3) = -P/tau * [ t*cos(delta), t*sin(delta), sin(sigma) ]                                   ! passive axis-angle pair so a minus sign in front
    ax(4) = alpha
    if (sign(1.0_pReal,alpha) < 0.0_pReal) ax = -ax                                                 ! ensure alpha is positive
  end if

end function eu2ax


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert axis-angle pair to unit quaternion.
!--------------------------------------------------------------------------------------------------
pure function ax2qu(ax) result(qu)

  real(pReal), intent(in), dimension(4) :: ax
  real(pReal),             dimension(4) :: qu

  real(pReal)                           :: c, s


  if (dEq0(ax(4))) then
    qu = [ 1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal ]
  else
    c = cos(ax(4)*0.5_pReal)
    s = sin(ax(4)*0.5_pReal)
    qu = [ c, ax(1)*s, ax(2)*s, ax(3)*s ]
  end if

end function ax2qu


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert axis-angle pair to orientation matrix.
!--------------------------------------------------------------------------------------------------
pure function ax2om(ax) result(om)

  real(pReal), intent(in), dimension(4)   :: ax
  real(pReal),             dimension(3,3) :: om

  real(pReal)                             :: q, c, s, omc


  c = cos(ax(4))
  s = sin(ax(4))
  omc = 1.0_pReal-c

  om(1,1) = ax(1)**2*omc + c
  om(2,2) = ax(2)**2*omc + c
  om(3,3) = ax(3)**2*omc + c

  q = omc*ax(1)*ax(2)
  om(1,2) = q + s*ax(3)
  om(2,1) = q - s*ax(3)

  q = omc*ax(2)*ax(3)
  om(2,3) = q + s*ax(1)
  om(3,2) = q - s*ax(1)

  q = omc*ax(3)*ax(1)
  om(3,1) = q + s*ax(2)
  om(1,3) = q - s*ax(2)

  if (P > 0.0_pReal) om = transpose(om)

end function ax2om


!--------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Convert axis-angle pair to Bunge Euler angles.
!--------------------------------------------------------------------------------------------------
pure function ax2eu(ax) result(eu)

  real(pReal), intent(in), dimension(4) :: ax
  real(pReal),             dimension(3) :: eu


  eu = om2eu(ax2om(ax))

end function ax2eu


!--------------------------------------------------------------------------------------------------
!> @brief Multiply two quaternions.
!--------------------------------------------------------------------------------------------------
pure function multiplyQuaternion(qu1,qu2)

  real(pReal), dimension(4), intent(in) :: qu1, qu2
  real(pReal), dimension(4) :: multiplyQuaternion


  multiplyQuaternion(1) = qu1(1)*qu2(1) - qu1(2)*qu2(2) -      qu1(3)*qu2(3) - qu1(4)*qu2(4)
  multiplyQuaternion(2) = qu1(1)*qu2(2) + qu1(2)*qu2(1) + P * (qu1(3)*qu2(4) - qu1(4)*qu2(3))
  multiplyQuaternion(3) = qu1(1)*qu2(3) + qu1(3)*qu2(1) + P * (qu1(4)*qu2(2) - qu1(2)*qu2(4))
  multiplyQuaternion(4) = qu1(1)*qu2(4) + qu1(4)*qu2(1) + P * (qu1(2)*qu2(3) - qu1(3)*qu2(2))

end function multiplyQuaternion


!--------------------------------------------------------------------------------------------------
!> @brief Calculate conjugate complex of a quaternion.
!--------------------------------------------------------------------------------------------------
pure function conjugateQuaternion(qu)

  real(pReal), dimension(4), intent(in) :: qu
  real(pReal), dimension(4) :: conjugateQuaternion


  conjugateQuaternion = [qu(1), -qu(2), -qu(3), -qu(4)]

end function conjugateQuaternion


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some rotations functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  type(tRotation)                 :: R
  real(pReal), dimension(4)       :: qu, ax
  real(pReal), dimension(3)       :: x, eu, v3
  real(pReal), dimension(3,3)     :: om, t33
  real(pReal), dimension(3,3,3,3) :: t3333
  real(pReal), dimension(6,6)     :: C
  real(pReal) :: A,B
  integer :: i


  do i = 1, 20

    if(i==1) then
      qu = [1.0_pReal, 0.0_pReal, 0.0_pReal, 0.0_pReal]
    elseif(i==2) then
      qu = [1.0_pReal,-0.0_pReal,-0.0_pReal,-0.0_pReal]
    elseif(i==3) then
      qu = [0.0_pReal, 1.0_pReal, 0.0_pReal, 0.0_pReal]
    elseif(i==4) then
      qu = [0.0_pReal,0.0_pReal,1.0_pReal,0.0_pReal]
    elseif(i==5) then
      qu = [0.0_pReal, 0.0_pReal, 0.0_pReal, 1.0_pReal]
    else
      call random_number(x)
      A = sqrt(x(3))
      B = sqrt(1-0_pReal -x(3))
      qu = [cos(TAU*x(1))*A,&
            sin(TAU*x(2))*B,&
            cos(TAU*x(2))*B,&
            sin(TAU*x(1))*A]
      if(qu(1)<0.0_pReal) qu = qu * (-1.0_pReal)
    end if


    if(.not. quaternion_equal(om2qu(qu2om(qu)),qu)) error stop 'om2qu2om'
    if(.not. quaternion_equal(eu2qu(qu2eu(qu)),qu)) error stop 'eu2qu2eu'
    if(.not. quaternion_equal(ax2qu(qu2ax(qu)),qu)) error stop 'ax2qu2ax'

    om = qu2om(qu)
    if(.not. quaternion_equal(om2qu(eu2om(om2eu(om))),qu)) error stop 'eu2om2eu'
    if(.not. quaternion_equal(om2qu(ax2om(om2ax(om))),qu)) error stop 'ax2om2ax'

    eu = qu2eu(qu)
    if(.not. quaternion_equal(eu2qu(ax2eu(eu2ax(eu))),qu)) error stop 'ax2eu2ax'

    call R%fromMatrix(om)

    call random_number(v3)
    if (any(dNeq(R%rotVector(R%rotVector(v3),active=.true.),v3,1.0e-12_pReal))) &
      error stop 'rotVector'

    call random_number(t33)
    if (any(dNeq(R%rotTensor2(R%rotTensor2(t33),active=.true.),t33,1.0e-12_pReal))) &
      error stop 'rotTensor2'

    call random_number(t3333)
    if (any(dNeq(R%rotTensor4(R%rotTensor4(t3333),active=.true.),t3333,1.0e-12_pReal))) &
      error stop 'rotTensor4'

    call random_number(C)
    C = C+transpose(C)
    if (any(dNeq(R%rotStiffness(C), &
                 math_3333toVoigt66_stiffness(R%rotate(math_Voigt66to3333_stiffness(C))),1.0e-12_pReal))) &
      error stop 'rotStiffness'

    call R%fromQuaternion(qu * (1.0_pReal + merge(+5.e-9_pReal,-5.e-9_pReal, mod(i,2) == 0)))       ! allow reasonable tolerance for ASCII/YAML

  end do

  contains

  pure recursive function quaternion_equal(qu1,qu2) result(ok)

    real(pReal), intent(in), dimension(4) :: qu1,qu2
    logical :: ok

    ok = all(dEq(qu1,qu2,1.0e-7_pReal))
    if(dEq0(qu1(1),1.0e-12_pReal)) &
      ok = ok .or. all(dEq(-1.0_pReal*qu1,qu2,1.0e-7_pReal))

  end function quaternion_equal

end subroutine selfTest


end module rotations
