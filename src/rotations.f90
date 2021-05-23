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

!---------------------------------------------------------------------------------------------------
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
!---------------------------------------------------------------------------------------------------

module rotations
  use IO
  use math

  implicit none
  private

  real(pReal), parameter :: P = -1.0_pReal                                                          !< parameter for orientation conversion.

  type, public :: rotation
    real(pReal), dimension(4) :: q
    contains
      procedure, public :: asQuaternion
      procedure, public :: asEulers
      procedure, public :: asAxisAngle
      procedure, public :: asRodrigues
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
      procedure, public  :: rotTensor4sym
      procedure, public  :: misorientation
      procedure, public  :: standardize
  end type rotation

  real(pReal), parameter :: &
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
    rotations_init, &
    eu2om

contains

!--------------------------------------------------------------------------------------------------
!> @brief Do self test.
!--------------------------------------------------------------------------------------------------
subroutine rotations_init

  print'(/,a)', ' <<<+-  rotations init  -+>>>'; flush(IO_STDOUT)

  print*, 'D. Rowenhorst et al., Modelling and Simulation in Materials Science and Engineering 23:083501, 2015'
  print*, 'https://doi.org/10.1088/0965-0393/23/8/083501'

  call selfTest

end subroutine rotations_init


!---------------------------------------------------------------------------------------------------
! Return rotation in different representations
!---------------------------------------------------------------------------------------------------
pure function asQuaternion(self)

  class(rotation), intent(in) :: self
  real(pReal), dimension(4)   :: asQuaternion

  asQuaternion = self%q

end function asQuaternion
!---------------------------------------------------------------------------------------------------
pure function asEulers(self)

  class(rotation), intent(in) :: self
  real(pReal), dimension(3)   :: asEulers

  asEulers = qu2eu(self%q)

end function asEulers
!---------------------------------------------------------------------------------------------------
pure function asAxisAngle(self)

  class(rotation), intent(in) :: self
  real(pReal), dimension(4)   :: asAxisAngle

  asAxisAngle = qu2ax(self%q)

end function asAxisAngle
!---------------------------------------------------------------------------------------------------
pure function asMatrix(self)

  class(rotation), intent(in) :: self
  real(pReal), dimension(3,3) :: asMatrix

  asMatrix = qu2om(self%q)

end function asMatrix
!---------------------------------------------------------------------------------------------------
pure function asRodrigues(self)

  class(rotation), intent(in) :: self
  real(pReal), dimension(4)   :: asRodrigues

  asRodrigues = qu2ro(self%q)

end function asRodrigues
!---------------------------------------------------------------------------------------------------
pure function asHomochoric(self)

  class(rotation), intent(in) :: self
  real(pReal), dimension(3)   :: asHomochoric

  asHomochoric = qu2ho(self%q)

end function asHomochoric

!---------------------------------------------------------------------------------------------------
! Initialize rotation from different representations
!---------------------------------------------------------------------------------------------------
subroutine fromQuaternion(self,qu)

  class(rotation), intent(out)          :: self
  real(pReal), dimension(4), intent(in) :: qu

  if (dNeq(norm2(qu),1.0_pReal,1.0e-8_pReal)) call IO_error(402,ext_msg='fromQuaternion')

  self%q = qu

end subroutine fromQuaternion
!---------------------------------------------------------------------------------------------------
subroutine fromEulers(self,eu,degrees)

  class(rotation), intent(out)          :: self
  real(pReal), dimension(3), intent(in) :: eu
  logical, intent(in), optional         :: degrees

  real(pReal), dimension(3)             :: Eulers

  if (.not. present(degrees)) then
    Eulers = eu
  else
    Eulers = merge(eu*INRAD,eu,degrees)
  endif

  if (any(Eulers<0.0_pReal) .or. any(Eulers>2.0_pReal*PI) .or. Eulers(2) > PI) &
    call IO_error(402,ext_msg='fromEulers')

  self%q = eu2qu(Eulers)

end subroutine fromEulers
!---------------------------------------------------------------------------------------------------
subroutine fromAxisAngle(self,ax,degrees,P)

  class(rotation), intent(out)          :: self
  real(pReal), dimension(4), intent(in) :: ax
  logical, intent(in), optional         :: degrees
  integer, intent(in), optional         :: P

  real(pReal)                           :: angle
  real(pReal),dimension(3)              :: axis

  if (.not. present(degrees)) then
    angle = ax(4)
  else
    angle = merge(ax(4)*INRAD,ax(4),degrees)
  endif

  if (.not. present(P)) then
    axis = ax(1:3)
  else
    axis = ax(1:3) * merge(-1.0_pReal,1.0_pReal,P == 1)
    if(abs(P) /= 1) call IO_error(402,ext_msg='fromAxisAngle (P)')
  endif

  if(dNeq(norm2(axis),1.0_pReal) .or. angle < 0.0_pReal .or. angle > PI) &
    call IO_error(402,ext_msg='fromAxisAngle')

  self%q = ax2qu([axis,angle])

end subroutine fromAxisAngle
!---------------------------------------------------------------------------------------------------
subroutine fromMatrix(self,om)

  class(rotation), intent(out)            :: self
  real(pReal), dimension(3,3), intent(in) :: om

  if (dNeq(math_det33(om),1.0_pReal,tol=1.0e-5_pReal)) &
    call IO_error(402,ext_msg='fromMatrix')

  self%q = om2qu(om)

end subroutine fromMatrix
!---------------------------------------------------------------------------------------------------


!---------------------------------------------------------------------------------------------------
!> @brief: Rotate a rotation
!---------------------------------------------------------------------------------------------------
pure elemental function rotRot__(self,R) result(rRot)

  type(rotation)              :: rRot
  class(rotation), intent(in) :: self,R

  rRot = rotation(multiply_quaternion(self%q,R%q))
  call rRot%standardize()

end function rotRot__


!---------------------------------------------------------------------------------------------------
!> @brief quaternion representation with positive q
!---------------------------------------------------------------------------------------------------
pure elemental subroutine standardize(self)

  class(rotation), intent(inout) :: self

  if (self%q(1) < 0.0_pReal) self%q = - self%q

end subroutine standardize


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Rotate a vector passively (default) or actively.
!---------------------------------------------------------------------------------------------------
pure function rotVector(self,v,active) result(vRot)

  real(pReal),                 dimension(3) :: vRot
  class(rotation), intent(in)               :: self
  real(pReal),     intent(in), dimension(3) :: v
  logical,         intent(in), optional     :: active

  real(pReal), dimension(4) :: v_normed, q
  logical                   :: passive

  if (present(active)) then
    passive = .not. active
  else
    passive = .true.
  endif

  if (dEq0(norm2(v))) then
    vRot = v
  else
    v_normed = [0.0_pReal,v]/norm2(v)
    if (passive) then
      q = multiply_quaternion(self%q, multiply_quaternion(v_normed, conjugate_quaternion(self%q)))
    else
      q = multiply_quaternion(conjugate_quaternion(self%q), multiply_quaternion(v_normed, self%q))
    endif
    vRot = q(2:4)*norm2(v)
  endif

end function rotVector


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Rotate a rank-2 tensor passively (default) or actively.
!> @details: Rotation is based on rotation matrix
!---------------------------------------------------------------------------------------------------
pure function rotTensor2(self,T,active) result(tRot)

  real(pReal),                 dimension(3,3) :: tRot
  class(rotation), intent(in)                 :: self
  real(pReal),     intent(in), dimension(3,3) :: T
  logical,         intent(in), optional       :: active

  logical           :: passive

  if (present(active)) then
    passive = .not. active
  else
    passive = .true.
  endif

  if (passive) then
    tRot = matmul(matmul(self%asMatrix(),T),transpose(self%asMatrix()))
  else
    tRot = matmul(matmul(transpose(self%asMatrix()),T),self%asMatrix())
  endif

end function rotTensor2


!---------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief rotate a rank-4 tensor passively (default) or actively
!> @details: rotation is based on rotation matrix
!! ToDo: Need to check active/passive !!!
!---------------------------------------------------------------------------------------------------
pure function rotTensor4(self,T,active) result(tRot)

  real(pReal),                 dimension(3,3,3,3) :: tRot
  class(rotation), intent(in)                     :: self
  real(pReal),     intent(in), dimension(3,3,3,3) :: T
  logical,         intent(in), optional           :: active

  real(pReal), dimension(3,3) :: R
  integer :: i,j,k,l,m,n,o,p

  if (present(active)) then
    R = merge(transpose(self%asMatrix()),self%asMatrix(),active)
  else
    R = self%asMatrix()
  endif

  tRot = 0.0_pReal
  do i = 1,3;do j = 1,3;do k = 1,3;do l = 1,3
  do m = 1,3;do n = 1,3;do o = 1,3;do p = 1,3
    tRot(i,j,k,l) = tRot(i,j,k,l) &
                  + R(i,m) * R(j,n) * R(k,o) * R(l,p) * T(m,n,o,p)
  enddo; enddo; enddo; enddo; enddo; enddo; enddo; enddo

end function rotTensor4


!---------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief rotate a symmetric rank-4 tensor stored as (6,6) passively (default) or actively
!! ToDo: Need to check active/passive !!!
!---------------------------------------------------------------------------------------------------
pure function rotTensor4sym(self,T,active) result(tRot)

  real(pReal),                 dimension(6,6) :: tRot
  class(rotation), intent(in)                 :: self
  real(pReal),     intent(in), dimension(6,6) :: T
  logical,         intent(in), optional       :: active

  if (present(active)) then
    tRot = math_sym3333to66(rotTensor4(self,math_66toSym3333(T),active))
  else
    tRot = math_sym3333to66(rotTensor4(self,math_66toSym3333(T)))
  endif

end function rotTensor4sym


!---------------------------------------------------------------------------------------------------
!> @brief misorientation
!---------------------------------------------------------------------------------------------------
pure elemental function misorientation(self,other)

  type(rotation)              :: misorientation
  class(rotation), intent(in) :: self, other

  misorientation%q = multiply_quaternion(other%q, conjugate_quaternion(self%q))

end function misorientation


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to rotation matrix
!---------------------------------------------------------------------------------------------------
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

  if (P < 0.0_pReal) om = transpose(om)

end function qu2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to Euler angles
!---------------------------------------------------------------------------------------------------
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
  endif degenerated
  where(eu<0.0_pReal) eu = mod(eu+2.0_pReal*PI,[2.0_pReal*PI,PI,2.0_pReal*PI])

end function qu2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to axis angle pair
!---------------------------------------------------------------------------------------------------
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


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function qu2ro(qu) result(ro)

  real(pReal), intent(in), dimension(4) :: qu
  real(pReal),             dimension(4) :: ro

  real(pReal)                           :: s
  real(pReal), parameter                :: thr = 1.0e-8_pReal

  if (abs(qu(1)) < thr) then
    ro =   [qu(2),  qu(3),  qu(4), IEEE_value(1.0_pReal,IEEE_positive_inf)]
  else
    s = norm2(qu(2:4))
    if (s < thr) then
      ro = [0.0_pReal, 0.0_pReal, P, 0.0_pReal]
    else
      ro = [qu(2)/s,qu(3)/s,qu(4)/s, tan(acos(math_clip(qu(1),-1.0_pReal,1.0_pReal)))]
    endif

  end if

end function qu2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to homochoric
!---------------------------------------------------------------------------------------------------
pure function qu2ho(qu) result(ho)

  real(pReal), intent(in), dimension(4) :: qu
  real(pReal),             dimension(3) :: ho

  real(pReal)                           :: omega, f

  omega = 2.0 * acos(math_clip(qu(1),-1.0_pReal,1.0_pReal))

  if (dEq0(omega,tol=1.e-5_pReal)) then
    ho = [ 0.0_pReal, 0.0_pReal, 0.0_pReal ]
  else
    ho = qu(2:4)
    f  = 0.75_pReal * ( omega - sin(omega) )
    ho = ho/norm2(ho)* f**(1.0_pReal/3.0_pReal)
  end if

end function qu2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert unit quaternion to cubochoric
!---------------------------------------------------------------------------------------------------
pure function qu2cu(qu) result(cu)

  real(pReal), intent(in), dimension(4) :: qu
  real(pReal),             dimension(3) :: cu

  cu = ho2cu(qu2ho(qu))

end function qu2cu


!---------------------------------------------------------------------------------------------------
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief convert rotation matrix to cubochoric
!> @details the original formulation (direct conversion) had (numerical?) issues
!---------------------------------------------------------------------------------------------------
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
      endif
  endif
  if(qu(1)<0._pReal) qu =-1.0_pReal * qu
  qu = qu*[1.0_pReal,P,P,P]

end function om2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief orientation matrix to Euler angles
!> @details Two step check for special cases to avoid invalid operations (not needed for python)
!---------------------------------------------------------------------------------------------------
pure function om2eu(om) result(eu)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(3)   :: eu
  real(pReal)                             :: zeta

  if    (dNeq(abs(om(3,3)),1.0_pReal,1.e-8_pReal)) then
    zeta = 1.0_pReal/sqrt(math_clip(1.0_pReal-om(3,3)**2.0_pReal,1e-64_pReal,1.0_pReal))
    eu = [atan2(om(3,1)*zeta,-om(3,2)*zeta), &
          acos(math_clip(om(3,3),-1.0_pReal,1.0_pReal)), &
          atan2(om(1,3)*zeta, om(2,3)*zeta)]
  else
    eu = [atan2(om(1,2),om(1,1)), 0.5_pReal*PI*(1.0_pReal-om(3,3)),0.0_pReal ]
  end if
  where(abs(eu) < 1.e-8_pReal) eu = 0.0_pReal
  where(eu<0.0_pReal)          eu = mod(eu+2.0_pReal*PI,[2.0_pReal*PI,PI,2.0_pReal*PI])

end function om2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert orientation matrix to axis angle pair
!---------------------------------------------------------------------------------------------------
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
#if defined(__GFORTRAN__) &&  __GNUC__<9
    i = maxloc(merge(1,0,cEq(cmplx(Wr,Wi,pReal),cmplx(1.0_pReal,0.0_pReal,pReal),tol=1.0e-14_pReal)),dim=1)
#else
    i = findloc(cEq(cmplx(Wr,Wi,pReal),cmplx(1.0_pReal,0.0_pReal,pReal),tol=1.0e-14_pReal),.true.,dim=1) !find eigenvalue (1,0)
#endif
    if (i == 0) error stop 'om2ax conversion failed'
    ax(1:3) = VR(1:3,i)
    where (                dNeq0([om(2,3)-om(3,2), om(3,1)-om(1,3), om(1,2)-om(2,1)])) &
      ax(1:3) = sign(ax(1:3),-P *[om(2,3)-om(3,2), om(3,1)-om(1,3), om(1,2)-om(2,1)])
  endif

end function om2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function om2ro(om) result(ro)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(4)   :: ro

  ro = eu2ro(om2eu(om))

end function om2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to homochoric
!---------------------------------------------------------------------------------------------------
function om2ho(om) result(ho)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(3)   :: ho

  ho = ax2ho(om2ax(om))

end function om2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert rotation matrix to cubochoric
!---------------------------------------------------------------------------------------------------
function om2cu(om) result(cu)

  real(pReal), intent(in), dimension(3,3) :: om
  real(pReal),             dimension(3)   :: cu

  cu = ho2cu(om2ho(om))

end function om2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Euler angles to unit quaternion
!---------------------------------------------------------------------------------------------------
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
  if(qu(1) < 0.0_pReal) qu = qu * (-1.0_pReal)

end function eu2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Euler angles to orientation matrix
!---------------------------------------------------------------------------------------------------
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


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert euler to axis angle
!---------------------------------------------------------------------------------------------------
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
    if (alpha < 0.0_pReal) ax = -ax                                                                 ! ensure alpha is positive
  end if

end function eu2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief Euler angles to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function eu2ro(eu) result(ro)

  real(pReal), intent(in), dimension(3) :: eu
  real(pReal),             dimension(4) :: ro

  ro = eu2ax(eu)
  if (ro(4) >= PI) then
    ro(4) = IEEE_value(ro(4),IEEE_positive_inf)
  elseif(dEq0(ro(4))) then
    ro = [ 0.0_pReal, 0.0_pReal, P, 0.0_pReal ]
  else
    ro(4) = tan(ro(4)*0.5_pReal)
  end if

end function eu2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Euler angles to homochoric
!---------------------------------------------------------------------------------------------------
pure function eu2ho(eu) result(ho)

  real(pReal), intent(in), dimension(3) :: eu
  real(pReal),             dimension(3) :: ho

  ho = ax2ho(eu2ax(eu))

end function eu2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Euler angles to cubochoric
!---------------------------------------------------------------------------------------------------
function eu2cu(eu) result(cu)

  real(pReal), intent(in), dimension(3) :: eu
  real(pReal),             dimension(3) :: cu

  cu = ho2cu(eu2ho(eu))

end function eu2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to quaternion
!---------------------------------------------------------------------------------------------------
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


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to orientation matrix
!---------------------------------------------------------------------------------------------------
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


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to Euler angles
!---------------------------------------------------------------------------------------------------
pure function ax2eu(ax) result(eu)

  real(pReal), intent(in), dimension(4) :: ax
  real(pReal),             dimension(3) :: eu

  eu = om2eu(ax2om(ax))

end function ax2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function ax2ro(ax) result(ro)

  real(pReal), intent(in), dimension(4) :: ax
  real(pReal),             dimension(4) :: ro

  real(pReal), parameter                :: thr = 1.0e-7_pReal

  if (dEq0(ax(4))) then
    ro = [ 0.0_pReal, 0.0_pReal, P, 0.0_pReal ]
  else
    ro(1:3) =  ax(1:3)
    ! we need to deal with the 180 degree case
    ro(4) = merge(IEEE_value(ro(4),IEEE_positive_inf),tan(ax(4)*0.5_pReal),abs(ax(4)-PI) < thr)
  end if

end function ax2ro


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to homochoric
!---------------------------------------------------------------------------------------------------
pure function ax2ho(ax) result(ho)

  real(pReal), intent(in), dimension(4) :: ax
  real(pReal),             dimension(3) :: ho

  real(pReal)                           :: f

  f = 0.75_pReal * ( ax(4) - sin(ax(4)) )
  f = f**(1.0_pReal/3.0_pReal)
  ho = ax(1:3) * f

end function ax2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert axis angle pair to cubochoric
!---------------------------------------------------------------------------------------------------
function ax2cu(ax) result(cu)

  real(pReal), intent(in), dimension(4) :: ax
  real(pReal),             dimension(3) :: cu

  cu = ho2cu(ax2ho(ax))

end function ax2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to unit quaternion
!---------------------------------------------------------------------------------------------------
pure function ro2qu(ro) result(qu)

  real(pReal), intent(in), dimension(4) :: ro
  real(pReal),             dimension(4) :: qu

  qu = ax2qu(ro2ax(ro))

end function ro2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to rotation matrix
!---------------------------------------------------------------------------------------------------
pure function ro2om(ro) result(om)

  real(pReal), intent(in), dimension(4)   :: ro
  real(pReal),             dimension(3,3) :: om

  om = ax2om(ro2ax(ro))

end function ro2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to Euler angles
!---------------------------------------------------------------------------------------------------
pure function ro2eu(ro) result(eu)

  real(pReal), intent(in), dimension(4) :: ro
  real(pReal),             dimension(3) :: eu

  eu = om2eu(ro2om(ro))

end function ro2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to axis angle pair
!---------------------------------------------------------------------------------------------------
pure function ro2ax(ro) result(ax)

  real(pReal), intent(in), dimension(4) :: ro
  real(pReal),             dimension(4) :: ax

  real(pReal)                           :: ta, angle

  ta = ro(4)

  if (.not. IEEE_is_finite(ta)) then
    ax = [ ro(1), ro(2), ro(3), PI ]
  elseif (dEq0(ta))  then
    ax = [ 0.0_pReal, 0.0_pReal, 1.0_pReal, 0.0_pReal ]
  else
    angle = 2.0_pReal*atan(ta)
    ta = 1.0_pReal/norm2(ro(1:3))
    ax = [ ro(1)/ta, ro(2)/ta, ro(3)/ta, angle ]
  end if

end function ro2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to homochoric
!---------------------------------------------------------------------------------------------------
pure function ro2ho(ro) result(ho)

  real(pReal), intent(in), dimension(4) :: ro
  real(pReal),             dimension(3) :: ho

  real(pReal)                           :: f

  if (dEq0(norm2(ro(1:3)))) then
    ho = [ 0.0_pReal, 0.0_pReal, 0.0_pReal ]
  else
    f = merge(2.0_pReal*atan(ro(4)) - sin(2.0_pReal*atan(ro(4))),PI, IEEE_is_finite(ro(4)))
    ho = ro(1:3) * (0.75_pReal*f)**(1.0_pReal/3.0_pReal)
  end if

end function ro2ho


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert Rodrigues vector to cubochoric
!---------------------------------------------------------------------------------------------------
pure function ro2cu(ro) result(cu)

  real(pReal), intent(in), dimension(4) :: ro
  real(pReal),             dimension(3) :: cu

  cu = ho2cu(ro2ho(ro))

end function ro2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to unit quaternion
!---------------------------------------------------------------------------------------------------
pure function ho2qu(ho) result(qu)

  real(pReal), intent(in), dimension(3) :: ho
  real(pReal),             dimension(4) :: qu

  qu = ax2qu(ho2ax(ho))

end function ho2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to rotation matrix
!---------------------------------------------------------------------------------------------------
pure function ho2om(ho) result(om)

  real(pReal), intent(in), dimension(3)   :: ho
  real(pReal),             dimension(3,3) :: om

  om = ax2om(ho2ax(ho))

end function ho2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to Euler angles
!---------------------------------------------------------------------------------------------------
pure function ho2eu(ho) result(eu)

  real(pReal), intent(in), dimension(3) :: ho
  real(pReal),             dimension(3) :: eu

  eu = ax2eu(ho2ax(ho))

end function ho2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert homochoric to axis angle pair
!---------------------------------------------------------------------------------------------------
pure function ho2ax(ho) result(ax)

  real(pReal), intent(in), dimension(3) :: ho
  real(pReal),             dimension(4) :: ax

  integer                               :: i
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
!> @brief convert homochoric to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function ho2ro(ho) result(ro)

  real(pReal), intent(in), dimension(3) :: ho
  real(pReal),             dimension(4) :: ro

  ro = ax2ro(ho2ax(ho))

end function ho2ro


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief convert homochoric to cubochoric
!--------------------------------------------------------------------------
pure function ho2cu(ho) result(cu)

  real(pReal), intent(in), dimension(3)   :: ho
  real(pReal),             dimension(3)   :: cu, xyz1, xyz3
  real(pReal),             dimension(2)   :: Tinv, xyz2
  real(pReal)                             :: rs, qxy, q2, sq2, q, tt
  integer,                 dimension(3,2) :: p

  rs = norm2(ho)
  if (rs > R1+1.e-6_pReal) then
    cu = IEEE_value(cu,IEEE_positive_inf)
    return
  endif

  center: if (all(dEq0(ho))) then
    cu = 0.0_pReal
  else center
    p = GetPyramidOrder(ho)
    xyz3 = ho(p(:,1))

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

    ! reverse the coordinates back to order according to the original pyramid number
    cu = xyz1(p(:,2))

  endif center

end function ho2cu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to unit quaternion
!---------------------------------------------------------------------------------------------------
pure function cu2qu(cu) result(qu)

  real(pReal), intent(in), dimension(3) :: cu
  real(pReal),             dimension(4) :: qu

  qu = ho2qu(cu2ho(cu))

end function cu2qu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to rotation matrix
!---------------------------------------------------------------------------------------------------
pure function cu2om(cu) result(om)

  real(pReal), intent(in), dimension(3)   :: cu
  real(pReal),             dimension(3,3) :: om

  om = ho2om(cu2ho(cu))

end function cu2om


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to Euler angles
!---------------------------------------------------------------------------------------------------
pure function cu2eu(cu) result(eu)

  real(pReal), intent(in), dimension(3) :: cu
  real(pReal),             dimension(3) :: eu

  eu = ho2eu(cu2ho(cu))

end function cu2eu


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to axis angle pair
!---------------------------------------------------------------------------------------------------
function cu2ax(cu) result(ax)

  real(pReal), intent(in), dimension(3) :: cu
  real(pReal),             dimension(4) :: ax

  ax = ho2ax(cu2ho(cu))

end function cu2ax


!---------------------------------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @brief convert cubochoric to Rodrigues vector
!---------------------------------------------------------------------------------------------------
pure function cu2ro(cu) result(ro)

  real(pReal), intent(in), dimension(3) :: cu
  real(pReal),             dimension(4) :: ro

  ro = ho2ro(cu2ho(cu))

end function cu2ro


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief map from 3D cubic grid to 3D ball
!--------------------------------------------------------------------------
pure function cu2ho(cu) result(ho)

  real(pReal), intent(in), dimension(3)   :: cu
  real(pReal),             dimension(3)   :: ho, LamXYZ, XYZ
  real(pReal),             dimension(2)   :: T
  real(pReal)                             :: c, s, q
  real(pReal), parameter                  :: eps = 1.0e-8_pReal
  integer,                 dimension(3,2) :: p
  integer,                 dimension(2)   :: order

  if (maxval(abs(cu)) > AP/2.0+eps) then
    ho = IEEE_value(cu,IEEE_positive_inf)
    return
  end if

  ! transform to the sphere grid via the curved square, and intercept the zero point
  center: if (all(dEq0(cu))) then
    ho  = 0.0_pReal
  else center
    ! get pyramide and scale by grid parameter ratio
    p = GetPyramidOrder(cu)
    XYZ = cu(p(:,1)) * sc

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

    ! reverse the coordinates back to order according to the original pyramid number
    ho = LamXYZ(p(:,2))

  endif center

end function cu2ho


!--------------------------------------------------------------------------
!> @author Marc De Graef, Carnegie Mellon University
!> @author Martin Diehl, Max-Planck-Institut für Eisenforschung GmbH
!> @brief Determine to which pyramid a point in a cubic grid belongs.
!--------------------------------------------------------------------------
pure function GetPyramidOrder(xyz)

  real(pReal),intent(in),dimension(3)   :: xyz
  integer,               dimension(3,2) :: GetPyramidOrder

  if      (((abs(xyz(1)) <=  xyz(3)).and.(abs(xyz(2)) <=  xyz(3))) .or. &
           ((abs(xyz(1)) <= -xyz(3)).and.(abs(xyz(2)) <= -xyz(3)))) then
    GetPyramidOrder = reshape([[1,2,3],[1,2,3]],[3,2])
  else if (((abs(xyz(3)) <=  xyz(1)).and.(abs(xyz(2)) <=  xyz(1))) .or. &
           ((abs(xyz(3)) <= -xyz(1)).and.(abs(xyz(2)) <= -xyz(1)))) then
    GetPyramidOrder = reshape([[2,3,1],[3,1,2]],[3,2])
  else if (((abs(xyz(1)) <=  xyz(2)).and.(abs(xyz(3)) <=  xyz(2))) .or. &
           ((abs(xyz(1)) <= -xyz(2)).and.(abs(xyz(3)) <= -xyz(2)))) then
    GetPyramidOrder = reshape([[3,1,2],[2,3,1]],[3,2])
  else
    GetPyramidOrder = -1                                                                            ! should be impossible, but might simplify debugging
  end if

end function GetPyramidOrder


!--------------------------------------------------------------------------------------------------
!> @brief Multiply two quaternions.
!--------------------------------------------------------------------------------------------------
pure function multiply_quaternion(qu1,qu2)

  real(pReal), dimension(4), intent(in) :: qu1, qu2
  real(pReal), dimension(4) :: multiply_quaternion


  multiply_quaternion(1) = qu1(1)*qu2(1) - qu1(2)*qu2(2) -      qu1(3)*qu2(3) - qu1(4)*qu2(4)
  multiply_quaternion(2) = qu1(1)*qu2(2) + qu1(2)*qu2(1) + P * (qu1(3)*qu2(4) - qu1(4)*qu2(3))
  multiply_quaternion(3) = qu1(1)*qu2(3) + qu1(3)*qu2(1) + P * (qu1(4)*qu2(2) - qu1(2)*qu2(4))
  multiply_quaternion(4) = qu1(1)*qu2(4) + qu1(4)*qu2(1) + P * (qu1(2)*qu2(3) - qu1(3)*qu2(2))

end function multiply_quaternion


!--------------------------------------------------------------------------------------------------
!> @brief Calculate conjugate complex of a quaternion.
!--------------------------------------------------------------------------------------------------
pure function conjugate_quaternion(qu)

  real(pReal), dimension(4), intent(in) :: qu
  real(pReal), dimension(4) :: conjugate_quaternion


  conjugate_quaternion = [qu(1), -qu(2), -qu(3), -qu(4)]


end function conjugate_quaternion


!--------------------------------------------------------------------------------------------------
!> @brief Check correctness of some rotations functions.
!--------------------------------------------------------------------------------------------------
subroutine selfTest()

  type(rotation)                  :: R
  real(pReal), dimension(4)       :: qu, ax, ro
  real(pReal), dimension(3)       :: x, eu, ho, v3
  real(pReal), dimension(3,3)     :: om, t33
  real(pReal), dimension(3,3,3,3) :: t3333
  real    :: A,B
  integer :: i


  do i = 1, 20

#if defined(__GFORTRAN__) && __GNUC__<9
    if(i<7) cycle
#endif

    if(i==1) then
      qu = om2qu(math_I3)
    elseif(i==2) then
      qu = eu2qu([0.0_pReal,0.0_pReal,0.0_pReal])
    elseif(i==3) then
      qu = eu2qu([2.0_pReal*PI,PI,2.0_pReal*PI])
    elseif(i==4) then
      qu = [0.0_pReal,0.0_pReal,1.0_pReal,0.0_pReal]
    elseif(i==5) then
      qu = ro2qu([1.0_pReal,0.0_pReal,0.0_pReal,IEEE_value(1.0_pReal, IEEE_positive_inf)])
    elseif(i==6) then
      qu = ax2qu([1.0_pReal,0.0_pReal,0.0_pReal,0.0_pReal])
    else
      call random_number(x)
      A = sqrt(x(3))
      B = sqrt(1-0_pReal -x(3))
      qu = [cos(2.0_pReal*PI*x(1))*A,&
            sin(2.0_pReal*PI*x(2))*B,&
            cos(2.0_pReal*PI*x(2))*B,&
            sin(2.0_pReal*PI*x(1))*A]
      if(qu(1)<0.0_pReal) qu = qu * (-1.0_pReal)
    endif

    if(.not. quaternion_equal(om2qu(qu2om(qu)),qu)) error stop 'om2qu/qu2om'
    if(.not. quaternion_equal(eu2qu(qu2eu(qu)),qu)) error stop 'eu2qu/qu2eu'
    if(.not. quaternion_equal(ax2qu(qu2ax(qu)),qu)) error stop 'ax2qu/qu2ax'
    if(.not. quaternion_equal(ro2qu(qu2ro(qu)),qu)) error stop 'ro2qu/qu2ro'
    if(.not. quaternion_equal(ho2qu(qu2ho(qu)),qu)) error stop 'ho2qu/qu2ho'
    if(.not. quaternion_equal(cu2qu(qu2cu(qu)),qu)) error stop 'cu2qu/qu2cu'

    om = qu2om(qu)
    if(.not. quaternion_equal(om2qu(eu2om(om2eu(om))),qu)) error stop 'eu2om/om2eu'
    if(.not. quaternion_equal(om2qu(ax2om(om2ax(om))),qu)) error stop 'ax2om/om2ax'
    if(.not. quaternion_equal(om2qu(ro2om(om2ro(om))),qu)) error stop 'ro2om/om2ro'
    if(.not. quaternion_equal(om2qu(ho2om(om2ho(om))),qu)) error stop 'ho2om/om2ho'
    if(.not. quaternion_equal(om2qu(cu2om(om2cu(om))),qu)) error stop 'cu2om/om2cu'

    eu = qu2eu(qu)
    if(.not. quaternion_equal(eu2qu(ax2eu(eu2ax(eu))),qu)) error stop 'ax2eu/eu2ax'
    if(.not. quaternion_equal(eu2qu(ro2eu(eu2ro(eu))),qu)) error stop 'ro2eu/eu2ro'
    if(.not. quaternion_equal(eu2qu(ho2eu(eu2ho(eu))),qu)) error stop 'ho2eu/eu2ho'
    if(.not. quaternion_equal(eu2qu(cu2eu(eu2cu(eu))),qu)) error stop 'cu2eu/eu2cu'

    ax = qu2ax(qu)
    if(.not. quaternion_equal(ax2qu(ro2ax(ax2ro(ax))),qu)) error stop 'ro2ax/ax2ro'
    if(.not. quaternion_equal(ax2qu(ho2ax(ax2ho(ax))),qu)) error stop 'ho2ax/ax2ho'
    if(.not. quaternion_equal(ax2qu(cu2ax(ax2cu(ax))),qu)) error stop 'cu2ax/ax2cu'

    ro = qu2ro(qu)
    if(.not. quaternion_equal(ro2qu(ho2ro(ro2ho(ro))),qu)) error stop 'ho2ro/ro2ho'
    if(.not. quaternion_equal(ro2qu(cu2ro(ro2cu(ro))),qu)) error stop 'cu2ro/ro2cu'

    ho = qu2ho(qu)
    if(.not. quaternion_equal(ho2qu(cu2ho(ho2cu(ho))),qu)) error stop 'cu2ho/ho2cu'

    call R%fromMatrix(om)

    call random_number(v3)
    if(all(dNeq(R%rotVector(R%rotVector(v3),active=.true.),v3,1.0e-12_pReal))) &
      error stop 'rotVector'

    call random_number(t33)
    if(all(dNeq(R%rotTensor2(R%rotTensor2(t33),active=.true.),t33,1.0e-12_pReal))) &
      error stop 'rotTensor2'

    call random_number(t3333)
    if(all(dNeq(R%rotTensor4(R%rotTensor4(t3333),active=.true.),t3333,1.0e-12_pReal))) &
      error stop 'rotTensor4'

    call R%fromQuaternion(qu * (1.0_pReal + merge(+5.e-9_pReal,-5.e-9_pReal, mod(i,2) == 0)))       ! allow reasonable tolerance for ASCII/YAML

  enddo

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
