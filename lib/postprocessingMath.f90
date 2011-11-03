!$Id$
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!all function below are taken from math.f90
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module math

real*8, parameter :: pi = 3.14159265358979323846264338327950288419716939937510

! *** 3x3 Identity ***
 real*8, dimension(3,3), parameter :: math_I3 = &
 reshape( (/ &
 1.0,0.0,0.0, &
 0.0,1.0,0.0, &
 0.0,0.0,1.0 /),(/3,3/))

contains
!**************************************************************************
! matrix multiplication 33x33 = 3x3
!**************************************************************************
pure function math_mul33x33(A,B)  

 implicit none

 integer  i,j
 real*8, dimension(3,3), intent(in) ::  A,B
 real*8, dimension(3,3) ::  math_mul33x33

 forall (i=1:3,j=1:3) math_mul33x33(i,j) = &
   A(i,1)*B(1,j) + A(i,2)*B(2,j) + A(i,3)*B(3,j)
 return

end function math_mul33x33

!**************************************************************************
! Cramer inversion of 3x3 matrix (subroutine)
!**************************************************************************
 PURE SUBROUTINE math_invert3x3(A, InvA, DetA, error)

!   Bestimmung der Determinanten und Inversen einer 3x3-Matrix
!   A      = Matrix A
!   InvA   = Inverse of A
!   DetA   = Determinant of A
!   error  = logical

 implicit none

 logical, intent(out) :: error

 real*8,dimension(3,3),intent(in)  :: A
 real*8,dimension(3,3),intent(out) :: InvA
 real*8, intent(out) :: DetA

 DetA =   A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) )&
        - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) )&
        + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) )

 if (DetA <= tiny(DetA)) then
   error = .true.
 else
   InvA(1,1) = (  A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / DetA
   InvA(2,1) = ( -A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / DetA
   InvA(3,1) = (  A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / DetA

   InvA(1,2) = ( -A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / DetA
   InvA(2,2) = (  A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / DetA
   InvA(3,2) = ( -A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / DetA

   InvA(1,3) = (  A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / DetA
   InvA(2,3) = ( -A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / DetA
   InvA(3,3) = (  A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / DetA
   
   error = .false.
 endif
 return

 END SUBROUTINE math_invert3x3
 
!********************************************************************
! determinant of a 3x3 matrix
!********************************************************************
 pure function math_det3x3(m)

 implicit none

 real*8, dimension(3,3), intent(in) :: m
 real*8 math_det3x3

 math_det3x3 = m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
              -m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
              +m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
 return

 end function math_det3x3
 
!****************************************************************
 pure subroutine math_pDecomposition(FE,U,R,error)
!-----FE = R.U 
!****************************************************************
 implicit none

 real*8, intent(in) :: FE(3,3)
 real*8, intent(out) :: R(3,3), U(3,3)
 logical, intent(out) :: error
 real*8 CE(3,3),EW1,EW2,EW3,EB1(3,3),EB2(3,3),EB3(3,3),UI(3,3),det

 error = .false.
 ce = math_mul33x33(transpose(FE),FE)

 CALL math_spectral1(CE,EW1,EW2,EW3,EB1,EB2,EB3)
 U=DSQRT(EW1)*EB1+DSQRT(EW2)*EB2+DSQRT(EW3)*EB3
 call math_invert3x3(U,UI,det,error)
 if (.not. error) R = math_mul33x33(FE,UI)

 return 
 
 end subroutine math_pDecomposition
 
!**************************************************************************
! Cramer inversion of 3x3 matrix (function)
!**************************************************************************
 pure function math_inv3x3(A)

!   direct Cramer inversion of matrix A.
!   returns all zeroes if not possible, i.e. if det close to zero

 implicit none

 real*8,dimension(3,3),intent(in)  :: A
 real*8 DetA

 real*8,dimension(3,3) :: math_inv3x3
 
 math_inv3x3 = 0.0

 DetA =   A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) )&
        - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) )&
        + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) )

 if (DetA > tiny(DetA)) then
   math_inv3x3(1,1) = (  A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / DetA
   math_inv3x3(2,1) = ( -A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / DetA
   math_inv3x3(3,1) = (  A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / DetA

   math_inv3x3(1,2) = ( -A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / DetA
   math_inv3x3(2,2) = (  A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / DetA
   math_inv3x3(3,2) = ( -A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / DetA

   math_inv3x3(1,3) = (  A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / DetA
   math_inv3x3(2,3) = ( -A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / DetA
   math_inv3x3(3,3) = (  A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / DetA
 endif
 return

 end function math_inv3x3
 
!********************************************************************** 
!  HAUPTINVARIANTEN HI1M, HI2M, HI3M DER 3X3 MATRIX M
!********************************************************************** 

 PURE SUBROUTINE math_hi(M,HI1M,HI2M,HI3M)
 implicit none

 real*8, intent(in) :: M(3,3) 
 real*8, intent(out) :: HI1M, HI2M, HI3M 

 HI1M=M(1,1)+M(2,2)+M(3,3)
 HI2M=HI1M**2/2.0-(M(1,1)**2+M(2,2)**2+M(3,3)**2)/2.0-M(1,2)*M(2,1)-M(1,3)*M(3,1)-M(2,3)*M(3,2)
 HI3M=math_det3x3(M)
! QUESTION: is 3rd equiv det(M) ?? if yes, use function math_det !agreed on YES
 return  

 END SUBROUTINE math_hi
 
!**********************************************************************
 pure subroutine math_spectral1(M,EW1,EW2,EW3,EB1,EB2,EB3)
!**** EIGENWERTE UND EIGENWERTBASIS DER SYMMETRISCHEN 3X3 MATRIX M

 implicit none

 real*8, intent(in) :: M(3,3)
 real*8, intent(out) :: EB1(3,3),EB2(3,3),EB3(3,3),EW1,EW2,EW3
 real*8 HI1M,HI2M,HI3M,TOL,R,S,T,P,Q,RHO,PHI,Y1,Y2,Y3,D1,D2,D3
 real*8 C1,C2,C3,M1(3,3),M2(3,3),M3(3,3),arg
 TOL=1.e-14
 CALL math_hi(M,HI1M,HI2M,HI3M)
 R=-HI1M
 S= HI2M
 T=-HI3M
 P=S-R**2.0/3.0
 Q=2.0/27.0*R**3.0-R*S/3.0+T
 EB1=0.0
 EB2=0.0
 EB3=0.0
 IF((ABS(P).LT.TOL).AND.(ABS(Q).LT.TOL))THEN
!   DREI GLEICHE EIGENWERTE
   EW1=HI1M/3.0
   EW2=EW1
   EW3=EW1
!   this is not really correct, but this way U is calculated
!   correctly in PDECOMPOSITION (correct is EB?=I)
   EB1(1,1)=1.0
   EB2(2,2)=1.0
   EB3(3,3)=1.0
 ELSE
   RHO=DSQRT(-3.0*P**3.0)/9.0
   arg=-Q/RHO/2.0
   if(arg.GT.1) arg=1
   if(arg.LT.-1) arg=-1
   PHI=DACOS(arg)
   Y1=2*RHO**(1.0/3.0)*DCOS(PHI/3.0)
   Y2=2*RHO**(1.0/3.0)*DCOS(PHI/3.0+2.0/3.0*PI)
   Y3=2*RHO**(1.0/3.0)*DCOS(PHI/3.0+4.0/3.0*PI)
   EW1=Y1-R/3.0
   EW2=Y2-R/3.0
   EW3=Y3-R/3.0
   C1=ABS(EW1-EW2)
   C2=ABS(EW2-EW3) 
   C3=ABS(EW3-EW1)

   IF(C1.LT.TOL) THEN
!  EW1 is equal to EW2
  D3=1.0/(EW3-EW1)/(EW3-EW2)
  M1=M-EW1*math_I3
  M2=M-EW2*math_I3
  EB3=math_mul33x33(M1,M2)*D3

  EB1=math_I3-EB3
!  both EB2 and EW2 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW2=0.0
   ELSE IF(C2.LT.TOL) THEN
!  EW2 is equal to EW3
  D1=1.0/(EW1-EW2)/(EW1-EW3)
  M2=M-math_I3*EW2
  M3=M-math_I3*EW3
  EB1=math_mul33x33(M2,M3)*D1
  EB2=math_I3-EB1
!  both EB3 and EW3 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW3=0.0
   ELSE IF(C3.LT.TOL) THEN
!  EW1 is equal to EW3
  D2=1.0/(EW2-EW1)/(EW2-EW3) 
  M1=M-math_I3*EW1
  M3=M-math_I3*EW3
  EB2=math_mul33x33(M1,M3)*D2
  EB1=math_I3-EB2
!  both EB3 and EW3 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW3=0.0
   ELSE
!  all three eigenvectors are different
  D1=1.0/(EW1-EW2)/(EW1-EW3)
  D2=1.0/(EW2-EW1)/(EW2-EW3) 
  D3=1.0/(EW3-EW1)/(EW3-EW2)
  M1=M-EW1*math_I3
  M2=M-EW2*math_I3
  M3=M-EW3*math_I3
  EB1=math_mul33x33(M2,M3)*D1
  EB2=math_mul33x33(M1,M3)*D2
  EB3=math_mul33x33(M1,M2)*D3

   END IF
 END IF
 RETURN
 END SUBROUTINE math_spectral1
 
!**************************************************************************
! volume of tetrahedron given by four vertices
!**************************************************************************
 pure function math_volTetrahedron(v1,v2,v3,v4)  

 implicit none

 real*8 math_volTetrahedron
 real*8, dimension (3), intent(in) :: v1,v2,v3,v4
 real*8, dimension (3,3) :: m

 m(:,1) = v1-v2
 m(:,2) = v2-v3
 m(:,3) = v3-v4

 math_volTetrahedron = math_det3x3(m)/6.0 

 end function math_volTetrahedron
 
!subroutines below are for postprocessing with python

!two small helper functions for indexing
! CAREFULL, index and location runs from 0 to N-1 (python style)

 function mesh_location(idx,resolution)
   integer, intent(in) :: idx
   integer, intent(in) :: resolution(3)
   integer :: mesh_location(3)
   mesh_location = (/modulo(idx/ resolution(3) / resolution(2),resolution(1)), &
                     modulo(idx/ resolution(3),                resolution(2)), &
                     modulo(idx,                               resolution(3))/)
 end function mesh_location
 
 function mesh_index(location,resolution)
   integer, intent(in) :: location(3)
   integer, intent(in) :: resolution(3)
   integer :: mesh_index
   
   mesh_index = modulo(location(3), resolution(3))     +&
               (modulo(location(2), resolution(2)))*resolution(3) +&
               (modulo(location(1), resolution(1)))*resolution(3)*resolution(2)
 end function mesh_index

 end module math
 

 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine mesh(res_x,res_y,res_z,geomdim,defgrad_av,centroids,nodes)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Routine to build a regular mesh of cubes for given coordinates (= center of the cubes)
!
 implicit none
 
 real*8 geomdim(3)
 integer res_x, res_y, res_z
 real*8 wrappedCentroids(res_x+2,res_y+2,res_z+2,3)
 real*8            nodes(res_x+1,res_y+1,res_z+1,3)
 real*8        centroids(res_x  ,res_y  ,res_z  ,3)

 integer, dimension(3,8) :: neighbor = reshape((/ &
                                     0, 0, 0, &
                                     1, 0, 0, &
                                     1, 1, 0, &
                                     0, 1, 0, &
                                     0, 0, 1, &
                                     1, 0, 1, &
                                     1, 1, 1, &
                                     0, 1, 1  &
                                    /), &
                                    (/3,8/))
 
 integer i,j,k,n
 real*8, dimension(3,3) :: defgrad_av
 integer, dimension(3) :: diag, shift, lookup, me, res
 
 nodes = 0.0
 diag = 1
 shift = 0
 lookup = 0
 
 res = (/res_x,res_y,res_z/)
 
 wrappedCentroids = 0.0
 wrappedCentroids(2:res_x+1,2:res_y+1,2:res_z+1,:) = centroids

 do k = 0,res_z+1
   do j = 0,res_y+1
     do i = 0,res_x+1
       if (k==0 .or. k==res_z+1 .or. &                               ! z skin
           j==0 .or. j==res_y+1 .or. &                               ! y skin
           i==0 .or. i==res_x+1      ) then                          ! x skin
         me = (/i,j,k/)                                              ! me on skin
         shift = sign(abs(res+diag-2*me)/(res+diag),res+diag-2*me)
         lookup = me-diag+shift*res
         wrappedCentroids(i+1,j+1,k+1,:) = centroids(lookup(1)+1,lookup(2)+1,lookup(3)+1,:) - &
                                           matmul(defgrad_av, shift*geomdim)
       endif
 enddo; enddo; enddo
 do k = 0,res_z
   do j = 0,res_y
     do i = 0,res_x
       do n = 1,8
         nodes(i+1,j+1,k+1,:) = nodes(i+1,j+1,k+1,:) + wrappedCentroids(i+1+neighbor(1,n), &
                                                                        j+1+neighbor(2,n), &
                                                                        k+1+neighbor(3,n), :)
 enddo; enddo; enddo; enddo
 nodes = nodes/8.0
  
end subroutine mesh
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine deformed(res_x,res_y,res_z,geomdim,defgrad,defgrad_av,coord_avgCorner)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Routine to calculate coordinates in current configuration for given defgrad
! using linear interpolation (blurres out high frequency defomation)
!
 implicit none
 real*8 geomdim(3)
 integer res_x, res_y, res_z
 real*8 coord(8,6,res_x,res_y,res_z,3)
 real*8 coord_avgOrder(8,res_x,res_y,res_z,3)
 real*8 coord_avgCorner(res_x,res_y,res_z,3)
 real*8 defgrad(res_x,res_y,res_z,3,3)
 integer, dimension(3,8) :: corner = reshape((/ &
                                              0, 0, 0,&
                                              1, 0, 0,&
                                              1, 1, 0,&
                                              0, 1, 0,&
                                              1, 1, 1,&
                                              0, 1, 1,&
                                              0, 0, 1,&
                                              1, 0, 1 &
                                                  /), &
                                               (/3,8/))
 integer, dimension(3,8) :: step = reshape((/ &
                                            1, 1, 1,&
                                           -1, 1, 1,& 
                                           -1,-1, 1,&
                                            1,-1, 1,&
                                           -1,-1,-1,&
                                            1,-1,-1,&
                                            1, 1,-1,&
                                           -1, 1,-1 & 
                                                /), &
                                             (/3,8/))
 integer, dimension(3,6) :: order = reshape((/ &
                                            1, 2, 3,&
                                            1, 3, 2,&
                                            2, 1, 3,&
                                            2, 3, 1,&
                                            3, 1, 2,&
                                            3, 2, 1 &
                                                /), &
                                             (/3,6/))

 real*8 myStep(3), fones(3), parameter_coords(3)
 real*8 defgrad_av(3,3)
 real*8 negative(3), positive(3)
 integer rear(3), init(3), ones(3), oppo(3), me(3), res(3)
 integer i, j, k, s, o 

 print*, 'Restore geometry using linear integration'
 print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
 print '(a,/,i5,i5,i5)', ' Resolution:', res_x,res_y,res_z
 ones = 1
 fones = 1.0
 coord_avgOrder=0.0
 
 res = (/res_x,res_y,res_z/)
 
 do s = 0, 7                               ! corners (from 0 to 7)
   init = corner(:,s+1)*(res-ones) +ones
   oppo = corner(:,mod((s+4),8)+1)*(res-ones) +ones
   do o=1,6 ! orders                      ! from 1 to 6)
     do k = init(order(3,o)), oppo(order(3,o)), step(order(3,o),s+1)
       rear(order(2,o)) = init(order(2,o))
       do j = init(order(2,o)), oppo(order(2,o)), step(order(2,o),s+1)
         rear(order(1,o)) = init(order(1,o))
         do i = init(order(1,o)), oppo(order(1,o)), step(order(1,o),s+1)
           me(order(1,o)) = i
           me(order(2,o)) = j
           me(order(3,o)) = k
           if ( (me(1)==init(1)).and.(me(2)==init(2)).and. (me(3)==init(3)) ) then
             coord(s+1,o,me(1),me(2),me(3),:) = geomdim * (matmul(defgrad_av,corner(:,s+1)) + &
                                                           matmul(defgrad(me(1),me(2),me(3),:,:),0.5*step(:,s+1)/res))

           else
             myStep = (me-rear)*geomdim/res
             coord(s+1,o,me(1),me(2),me(3),:) = coord(s+1,o,rear(1),rear(2),rear(3),:) + &
                                             0.5*matmul(defgrad(me(1),me(2),me(3),:,:) + &
                                                        defgrad(rear(1),rear(2),rear(3),:,:),myStep)
           endif
           rear = me
   enddo; enddo; enddo; enddo
   do i=1,6
     coord_avgOrder(s+1,:,:,:,:) = coord_avgOrder(s+1,:,:,:,:) + coord(s+1,i,:,:,:,:)/6.0
   enddo
 enddo

 do k=0, res_z-1
   do j=0, res_y-1
     do i=0, res_x-1
       parameter_coords = (2.0*(/i+0.0,j+0.0,k+0.0/)-real(res)+fones)/(real(res)-fones)
       positive = fones + parameter_coords
       negative = fones - parameter_coords
       coord_avgCorner(i+1,j+1,k+1,:) = ( coord_avgOrder(1,i+1,j+1,k+1,:) *negative(1)*negative(2)*negative(3)&
                                        + coord_avgOrder(2,i+1,j+1,k+1,:) *positive(1)*negative(2)*negative(3)&
                                        + coord_avgOrder(3,i+1,j+1,k+1,:) *positive(1)*positive(2)*negative(3)&
                                        + coord_avgOrder(4,i+1,j+1,k+1,:) *negative(1)*positive(2)*negative(3)&
                                        + coord_avgOrder(5,i+1,j+1,k+1,:) *positive(1)*positive(2)*positive(3)&
                                        + coord_avgOrder(6,i+1,j+1,k+1,:) *negative(1)*positive(2)*positive(3)&
                                        + coord_avgOrder(7,i+1,j+1,k+1,:) *negative(1)*negative(2)*positive(3)&
                                        + coord_avgOrder(8,i+1,j+1,k+1,:) *positive(1)*negative(2)*positive(3))*0.125
 enddo; enddo; enddo
end subroutine deformed

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine deformed_fft(res_x,res_y,res_z,geomdim,defgrad,defgrad_av,scaling,coords)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Routine to calculate coordinates in current configuration for given defgrad
! using integration in Fourier space (more accurate than deformed(...))
!
 implicit none
 integer res_x, res_y, res_z
 real*8 geomdim(3)
 real*8 defgrad(res_x,res_y,res_z,3,3)
 real*8 defgrad_av(3,3)
 real*8 scaling
 real*8 coords(res_x,res_y,res_z,3)
 complex*16 coords_fft(res_x/2+1,res_y,res_z,3)
 complex*16 defgrad_fft(res_x,res_y,res_z,3,3)
 integer i, j, k
 integer k_s(3)
 real*8 step(3)
 real*8 offset_coords(3)
 real*8, parameter :: pi = 3.14159265358979323846264338327950288419716939937510
 integer*8 ::  plan_fft(2)
 include 'fftw3.f' !header file for fftw3 (declaring variables). Library files are also needed
 
 print*, 'Restore geometry using FFT-based integration'
 print '(a,/,e12.5,e12.5,e12.5)', ' Dimension:', geomdim
 print '(a,/,i5,i5,i5)', ' Resolution:', res_x,res_y,res_z
 
 call dfftw_plan_many_dft(plan_fft(1),3,(/res_x,res_y,res_z/),9,&
   defgrad_fft,(/res_x,res_y,res_z/),1,res_x*res_y*res_z,&
   defgrad_fft,(/res_x,res_y,res_z/),1,res_x*res_y*res_z,FFTW_FORWARD,FFTW_PATIENT)  
    
 call dfftw_plan_many_dft_c2r(plan_fft(2),3,(/res_x,res_y,res_z/),3,&
   coords_fft,(/res_x/2+1,res_y,res_z/),1,(res_x/2+1)*res_y*res_z,&
   coords,    (/res_x,    res_y,res_z/),1, res_x*     res_y*res_z,FFTW_PATIENT) 
   
 coords_fft = 0.0
 defgrad_fft = defgrad

 step(1) = geomdim(1)/real(res_x)
 step(2) = geomdim(2)/real(res_y)
 step(3) = geomdim(3)/real(res_z)

 call dfftw_execute_dft(plan_fft(1), defgrad_fft, defgrad_fft)
 
 do k = 1, res_z
   k_s(3) = k-1
   if(k > res_z/2+1) k_s(3) = k_s(3)-res_z
   do j = 1, res_y
     k_s(2) = j-1
     if(j > res_y/2+1) k_s(2) = k_s(2)-res_y
     do i = 1, res_x/2+1
       k_s(1) = i-1
       if(i/=1) coords_fft(i,j,k,:) = coords_fft(i,j,k,:)&
                                    + defgrad_fft(i,j,k,:,1)*geomdim(1)/(real(k_s(1))*cmplx(0.0,1.0)*pi*2.0)
       if(j/=1) coords_fft(i,j,k,:) = coords_fft(i,j,k,:)&
                                    + defgrad_fft(i,j,k,:,2)*geomdim(2)/(real(k_s(2))*cmplx(0.0,1.0)*pi*2.0)
       if(k/=1) coords_fft(i,j,k,:) = coords_fft(i,j,k,:)&
                                    + defgrad_fft(i,j,k,:,3)*geomdim(3)/(real(k_s(3))*cmplx(0.0,1.0)*pi*2.0)
 enddo; enddo; enddo
 
 call dfftw_execute_dft_c2r(plan_fft(2), coords_fft, coords)
 coords = coords/real(res_x*res_y*res_z)

 offset_coords = matmul(defgrad(1,1,1,:,:),step/2.0) - scaling*coords(1,1,1,:)
 do k = 1, res_z; do j = 1, res_y; do i = 1, res_x
     coords(i,j,k,:) =  scaling*coords(i,j,k,:) + offset_coords + matmul(defgrad_av,&
                                                  (/step(1)*real(i-1),&
                                                    step(2)*real(j-1),&
                                                    step(3)*real(k-1)/))

 enddo; enddo; enddo
end subroutine deformed_fft


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine volume_compare(res_x,res_y,res_z,geomdim,nodes,defgrad,volume_mismatch)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Routine to calculate the mismatch between volume of reconstructed (compatible)
! cube and determinant of defgrad at the FP

 use math
 implicit none
 
 real*8 geomdim(3)
 integer res_x, res_y, res_z
 real*8            nodes(res_x+1,res_y+1,res_z+1,3)
 real*8          defgrad(res_x  ,res_y  ,res_z  ,3,3)
 real*8  volume_mismatch(res_x  ,res_y  ,res_z  )
 real*8  coords(8,3)
 integer i,j,k
 real*8 vol_initial

 print*, 'Calculating volume mismatch'
 vol_initial = geomdim(1)*geomdim(2)*geomdim(3)/real(res_x)/real(res_y)/real(res_z)
 do k = 1,res_z
   do j = 1,res_y
     do i = 1,res_x 
       coords(1,:) = nodes(i  ,j  ,k  ,:)
       coords(2,:) = nodes(i+1,j  ,k  ,:)
       coords(3,:) = nodes(i+1,j+1,k  ,:)
       coords(4,:) = nodes(i  ,j+1,k  ,:)
       coords(5,:) = nodes(i  ,j,  k+1,:)
       coords(6,:) = nodes(i+1,j  ,k+1,:)
       coords(7,:) = nodes(i+1,j+1,k+1,:)
       coords(8,:) = nodes(i  ,j+1,k+1,:)
       volume_mismatch(i,j,k) = abs(math_volTetrahedron(coords(7,:),coords(1,:),coords(8,:),coords(4,:))) &
                              + abs(math_volTetrahedron(coords(7,:),coords(1,:),coords(8,:),coords(5,:))) &
                              + abs(math_volTetrahedron(coords(7,:),coords(1,:),coords(3,:),coords(4,:))) &
                              + abs(math_volTetrahedron(coords(7,:),coords(1,:),coords(3,:),coords(2,:))) &
                              + abs(math_volTetrahedron(coords(7,:),coords(5,:),coords(2,:),coords(6,:))) &
                              + abs(math_volTetrahedron(coords(7,:),coords(5,:),coords(2,:),coords(1,:)))
       volume_mismatch(i,j,k) = volume_mismatch(i,j,k)/math_det3x3(defgrad(i,j,k,:,:))
 enddo; enddo; enddo
 volume_mismatch = volume_mismatch/vol_initial
end subroutine volume_compare

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine shape_compare(res_x,res_y,res_z,geomdim,nodes,centroids,defgrad,shape_mismatch)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
! Routine to calculate the mismatch between the vectors from the central point to
! the corners of reconstructed (combatible) volume element and the vectors calculated by deforming 
! the initial volume element with the  current deformation gradient 
 implicit none
 
 real*8 geomdim(3)
 integer res_x, res_y, res_z
 real*8            nodes(res_x+1,res_y+1,res_z+1,3)
 real*8        centroids(res_x  ,res_y  ,res_z  ,3)
 real*8          defgrad(res_x  ,res_y  ,res_z  ,3,3)
 real*8   shape_mismatch(res_x  ,res_y  ,res_z)
 real*8  coords_initial(8,3)
 integer i,j,k

 print*, 'Calculating shape mismatch'
 coords_initial(1,:) = (/-geomdim(1)/2.0/real(res_x),-geomdim(2)/2.0/real(res_y),-geomdim(3)/2.0/real(res_z)/)
 coords_initial(2,:) = (/+geomdim(1)/2.0/real(res_x),-geomdim(2)/2.0/real(res_y),-geomdim(3)/2.0/real(res_z)/)
 coords_initial(3,:) = (/+geomdim(1)/2.0/real(res_x),+geomdim(2)/2.0/real(res_y),-geomdim(3)/2.0/real(res_z)/)
 coords_initial(4,:) = (/-geomdim(1)/2.0/real(res_x),+geomdim(2)/2.0/real(res_y),-geomdim(3)/2.0/real(res_z)/)
 coords_initial(5,:) = (/-geomdim(1)/2.0/real(res_x),-geomdim(2)/2.0/real(res_y),+geomdim(3)/2.0/real(res_z)/)
 coords_initial(6,:) = (/+geomdim(1)/2.0/real(res_x),-geomdim(2)/2.0/real(res_y),+geomdim(3)/2.0/real(res_z)/)
 coords_initial(7,:) = (/+geomdim(1)/2.0/real(res_x),+geomdim(2)/2.0/real(res_y),+geomdim(3)/2.0/real(res_z)/)
 coords_initial(8,:) = (/-geomdim(1)/2.0/real(res_x),+geomdim(2)/2.0/real(res_y),+geomdim(3)/2.0/real(res_z)/)
 do i=1,8
 enddo
 do k = 1,res_z
   do j = 1,res_y
     do i = 1,res_x
       shape_mismatch(i,j,k) = &
           sqrt(sum((nodes(i  ,j  ,k  ,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(1,:)))**2.0))&
         + sqrt(sum((nodes(i+1,j  ,k  ,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(2,:)))**2.0))&
         + sqrt(sum((nodes(i+1,j+1,k  ,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(3,:)))**2.0))&
         + sqrt(sum((nodes(i  ,j+1,k  ,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(4,:)))**2.0))&
         + sqrt(sum((nodes(i  ,j,  k+1,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(5,:)))**2.0))&
         + sqrt(sum((nodes(i+1,j  ,k+1,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(6,:)))**2.0))&
         + sqrt(sum((nodes(i+1,j+1,k+1,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(7,:)))**2.0))&
         + sqrt(sum((nodes(i  ,j+1,k+1,:) - centroids(i,j,k,:) - matmul(defgrad(i,j,k,:,:), coords_initial(8,:)))**2.0))
 enddo; enddo; enddo
 end subroutine shape_compare

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine inverse_reconstruction(res_x,res_y,res_z,reference_configuration,current_configuration,defgrad)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Routine to calculate deformation gradient from reference and current configuration
! NOT WORKING BY NOW!!!!!!!!!!!!!
!
 use math
 implicit none
 integer res_x, res_y, res_z
 real*8 reference_configuration(res_x+1,res_y+1,res_z+1,3)
 real*8 current_configuration(res_x+1,res_y+1,res_z+1,3)
 real*8 defgrad(res_x,res_y,res_z,3,3)
 real*8 delta, tolerance, res, res_center
 real*8 reference(8,3)
 real*8 current(8,3)
 real*8 defgrad_temp(3,3)
 real*8 dres_dF(3,3)
 real*8 identity(3,3)
 real*8 ref_bar(3)
 real*8 current_bar(3)
 real*8 r(8)
 real*8 differentiate(9,3,3)
 integer i, j, k, m, l, x, y, o 
 
 identity = 0.0
 identity(1,1) = 1.0
 identity(2,2) = 1.0
 identity(3,3) = 1.0
 
 differentiate = 0.0
 
 tolerance = 1e-10
 delta = 1e-9
 
 k = 0
 do j = 1, 3; do i = 1, 3
   k = k+1
   differentiate(k,i,j) = 1.0
 enddo; enddo
 
 do k = 1, res_z
   do j = 1, res_y
     do i = 1, res_x
       reference(1,:) = reference_configuration(i  ,j  ,k  ,:)
       reference(2,:) = reference_configuration(i+1,j  ,k  ,:)
       reference(3,:) = reference_configuration(i+1,j+1,k  ,:)
       reference(4,:) = reference_configuration(i  ,j+1,k  ,:)
       reference(5,:) = reference_configuration(i  ,j  ,k+1,:)
       reference(6,:) = reference_configuration(i+1,j  ,k+1,:)
       reference(7,:) = reference_configuration(i+1,j+1,k+1,:)
       reference(8,:) = reference_configuration(i  ,j+1,k+1,:)
       current(1,:) = current_configuration(i  ,j  ,k  ,:)
       current(2,:) = current_configuration(i+1,j  ,k  ,:)
       current(3,:) = current_configuration(i+1,j+1,k  ,:)
       current(4,:) = current_configuration(i  ,j+1,k  ,:)
       current(5,:) = current_configuration(i  ,j  ,k+1,:)
       current(6,:) = current_configuration(i+1,j  ,k+1,:)
       current(7,:) = current_configuration(i+1,j+1,k+1,:)
       current(8,:) = current_configuration(i  ,j+1,k+1,:)
       
       do o=1,3
         ref_bar(o) = sum(reference(:,o))/8.0
         current_bar(o) = sum(current(:,o))/8.0
       enddo
       
       do o=1,8
         reference(o,:) = reference(o,:) -ref_bar
         current(o,:) = current(o,:) -current_bar
       enddo
       
       defgrad_temp = identity
       res_center = 2.0*tolerance
       o=0
       do while(res_center >= tolerance)
         o = o + 1
         do l = 1,8                                                             ! loop over corners
           r(l) = sqrt(sum((current(l,:)-matmul(defgrad_temp,reference(l,:)))**2))    ! corner distance
         enddo
         res_center = sum(r*r)                                                 ! current residuum
         print*, 'res_center', res_center
         m=0
         do y=1,3; do x=1,3                                                     ! numerical differentiation
           m = m+1
           do l = 1,8
             r(l) = sqrt(sum((current(l,:)-matmul((defgrad_temp+differentiate(m,:,:)*delta),reference(l,:)))**2))  ! corner distance
           enddo
           res = sum(r*r)
           print*,'res step', m, res
           dres_dF(x,y)  = (res-res_center)/delta
         enddo; enddo
         print*, 'dres_dF', dres_dF
         print*, 'deltadef', math_inv3x3(dres_dF)*res_center
         defgrad_temp = defgrad_temp - math_inv3x3(dres_dF)*res_center          ! Newton--Raphson
         print*, o, res_center
!         pause
       enddo
       defgrad(i,j,k,:,:) = defgrad_temp       
 enddo; enddo; enddo

end subroutine inverse_reconstruction

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine tensor_avg(res_x,res_y,res_z,tensor,avg)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate average of tensor field
!

 implicit none
 integer res_x, res_y, res_z
 real*8 tensor(res_x,res_y,res_z,3,3)
 real*8 avg(3,3)
 real*8 wgt
 integer m,n
 
 wgt = 1/real(res_x*res_y*res_z)

 do m = 1,3; do n = 1,3
    avg(m,n) = sum(tensor(:,:,:,m,n)) * wgt   
 enddo; enddo
end subroutine tensor_avg
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine logstrain_spat(res_x,res_y,res_z,defgrad,logstrain_field)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate logarithmic strain in spatial configuration for given defgrad field
!
 use math
 implicit none
 integer res_x, res_y, res_z
 integer i, j, k
 real*8 defgrad(res_x,res_y,res_z,3,3)
 real*8 logstrain_field(res_x,res_y,res_z,3,3)
 real*8 temp33_Real(3,3), temp33_Real2(3,3)
 real*8 eigenvectorbasis(3,3,3)
 real*8 eigenvalue(3)
 logical errmatinv
 
 do k = 1, res_z; do j = 1, res_y; do i = 1, res_x
   call math_pDecomposition(defgrad(i,j,k,:,:),temp33_Real2,temp33_Real,errmatinv)  !store R in temp33_Real
   temp33_Real2 = math_inv3x3(temp33_Real)
   temp33_Real = math_mul33x33(defgrad(i,j,k,:,:),temp33_Real2)       ! v = F o inv(R), store in temp33_Real2
   call math_spectral1(temp33_Real, eigenvalue(1),          eigenvalue(2),            eigenvalue(3),&
                               eigenvectorbasis(1,:,:), eigenvectorbasis(2,:,:), eigenvectorbasis(3,:,:))
   eigenvalue = log(sqrt(eigenvalue))
   logstrain_field(i,j,k,:,:) = eigenvalue(1)*eigenvectorbasis(1,:,:)+&
                                eigenvalue(2)*eigenvectorbasis(2,:,:)+&
                                eigenvalue(3)*eigenvectorbasis(3,:,:)
 enddo; enddo; enddo
 end subroutine logstrain_spat
 
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine logstrain_mat(res_x,res_y,res_z,defgrad,logstrain_field)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate logarithmic strain in material configuration for given defgrad field
!
 use math
 implicit none
 integer res_x, res_y, res_z
 integer i, j, k
 real*8 defgrad(res_x,res_y,res_z,3,3)
 real*8 logstrain_field(res_x,res_y,res_z,3,3)
 real*8 temp33_Real(3,3), temp33_Real2(3,3)
 real*8 eigenvectorbasis(3,3,3)
 real*8 eigenvalue(3)
 logical errmatinv
 
 do k = 1, res_z; do j = 1, res_y; do i = 1, res_x
   call math_pDecomposition(defgrad(i,j,k,:,:),temp33_Real,temp33_Real2,errmatinv)  !store U in temp33_Real
   call math_spectral1(temp33_Real, eigenvalue(1),          eigenvalue(2),            eigenvalue(3),&
                               eigenvectorbasis(1,:,:), eigenvectorbasis(2,:,:), eigenvectorbasis(3,:,:))
   eigenvalue = log(sqrt(eigenvalue))
   logstrain_field(i,j,k,:,:) = eigenvalue(1)*eigenvectorbasis(1,:,:)+&
                                eigenvalue(2)*eigenvectorbasis(2,:,:)+&
                                eigenvalue(3)*eigenvectorbasis(3,:,:)
 enddo; enddo; enddo
 end subroutine logstrain_mat
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine calculate_cauchy(res_x,res_y,res_z,defgrad,p_stress,c_stress)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate cauchy stress for given PK1 stress and defgrad field
!
 use math
 implicit none
 integer res_x, res_y, res_z
 integer i, j, k
 real*8 defgrad(res_x,res_y,res_z,3,3)
 real*8 p_stress(res_x,res_y,res_z,3,3)
 real*8 c_stress(res_x,res_y,res_z,3,3)
 real*8 jacobi
 c_stress = 0.0
 do k = 1, res_z; do j = 1, res_y; do i = 1, res_x
   jacobi = math_det3x3(defgrad(i,j,k,:,:))
   c_stress(i,j,k,:,:) = matmul(p_stress(i,j,k,:,:),transpose(defgrad(i,j,k,:,:)))/jacobi
 enddo; enddo; enddo
end subroutine calculate_cauchy

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
subroutine calculate_mises(res_x,res_y,res_z,tensor,vm)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!calculate von Mises equivalent of tensor field
!
 implicit none
 integer res_x, res_y, res_z
 integer i, j, k
 real*8 tensor(res_x,res_y,res_z,3,3)
 real*8 vm(res_x,res_y,res_z,1)
 real*8 deviator(3,3)
 real*8 delta(3,3)
 real*8 J_2
 
 delta =0.0
 delta(1,1) = 1.0
 delta(2,2) = 1.0
 delta(3,3) = 1.0
 do k = 1, res_z; do j = 1, res_y; do i = 1, res_x
   deviator = tensor(i,j,k,:,:) - 1.0/3.0*tensor(i,j,k,1,1)*tensor(i,j,k,2,2)*tensor(i,j,k,3,3)*delta
   J_2 = deviator(1,1)*deviator(2,2)&
       + deviator(2,2)*deviator(3,3)&
       + deviator(1,1)*deviator(3,3)&
       - (deviator(1,2))**2&
       - (deviator(2,3))**2&
       - (deviator(1,3))**2
   vm(i,j,k,:) = sqrt(3*J_2)
 enddo; enddo; enddo
end subroutine calculate_mises

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine divergence_fft(res_x,res_y,res_z,vec_tens,geomdim,field,divergence_field)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculates divergence field using integration in Fourier space
!use vec_tens to decide if tensor (3) or vector (1)

 implicit none
 integer res_x, res_y, res_z, vec_tens
 real*8 geomdim(3)
 real*8 field(res_x,res_y,res_z,vec_tens,3)
 real*8 field_copy(res_x,res_y,res_z,vec_tens,3)
 real*8 xi(res_x,res_y,res_z,3)
 real*8 divergence_field(res_x,res_y,res_z,vec_tens)
 complex*16 divergence_field_fft(res_x/2+1,res_y,res_z,vec_tens)
 complex*16 field_fft(res_x,res_y,res_z,vec_tens,3)
 complex*16 img
 integer i, j, k
 real*8, parameter :: pi = 3.14159265358979323846264338327950288419716939937510
 integer*8 ::  plan_fft(2)
 include 'fftw3.f' !header file for fftw3 (declaring variables). Library files are also needed

 img = cmplx(0.0,1.0)

 call dfftw_plan_many_dft_r2c(plan_fft(1),3,(/res_x,res_y,res_z/),vec_tens*3,&
    field_copy,(/res_x,res_y,res_z/),1,res_x*res_y*res_z,&
    field_fft,(/res_x/2+1,res_y,res_z/),1,(res_x/2+1)*res_y*res_z,FFTW_PATIENT)

 call dfftw_plan_many_dft_c2r(plan_fft(2),3,(/res_x,res_y,res_z/),vec_tens,&
   divergence_field_fft,(/res_x/2+1,res_y,res_z/),1,(res_x/2+1)*res_y*res_z,&
   divergence_field,(/res_x,res_y,res_z/),1,res_x*res_y*res_z,FFTW_PATIENT)

! field_copy is destroyed during plan creation
 field_copy = field

 call dfftw_execute_dft_r2c(plan_fft(1), field_copy, field_fft)

 xi = 0.0
! Alternative calculation of discrete frequencies k_s, ordered as in FFTW (wrap around) 
! do k = 0,res_z/2 -1                                       
  ! do j = 0,res_y/2 -1
     ! do i = 0,res_x/2 -1
         ! xi(1+mod(res_x-i,res_x),1+mod(res_y-j,res_y),1+mod(res_z-k,res_z),:) = (/-i,-j,-k/)/geomdim
         ! xi(1+i,                 1+mod(res_y-j,res_y),1+mod(res_z-k,res_z),:) = (/ i,-j,-k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+j,                 1+mod(res_z-k,res_z),:) = (/-i, j,-k/)/geomdim
         ! xi(1+i,                 1+j,                 1+mod(res_z-k,res_z),:) = (/ i, j,-k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+mod(res_y-j,res_y),1+k,                 :) = (/-i,-j, k/)/geomdim
         ! xi(1+i,                 1+mod(res_y-j,res_y),1+k,                 :) = (/ i,-j, k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+j,                 1+k,                 :) = (/-i, j, k/)/geomdim
         ! xi(1+i,                 1+j,                 1+k,                 :) = (/ i, j, k/)/geomdim
         ! xi(1+i,                 1+j,                 1+k,                 :) = (/ i, j, k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+j,                 1+k,                 :) = (/-i, j, k/)/geomdim
         ! xi(1+i,                 1+mod(res_y-j,res_y),1+k,                 :) = (/ i,-j, k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+mod(res_y-j,res_y),1+k,                 :) = (/-i,-j, k/)/geomdim
         ! xi(1+i,                 1+j,                 1+mod(res_z-k,res_z),:) = (/ i, j,-k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+j,                 1+mod(res_z-k,res_z),:) = (/-i, j,-k/)/geomdim
         ! xi(1+i,                 1+mod(res_y-j,res_y),1+mod(res_z-k,res_z),:) = (/ i,-j,-k/)/geomdim
         ! xi(1+mod(res_x-i,res_x),1+mod(res_y-j,res_y),1+mod(res_z-k,res_z),:) = (/-i,-j,-k/)/geomdim
  ! enddo; enddo; enddo
  
 do k = 0, res_z-1         
   do j = 0, res_y-1
     do i = 0, res_x/2
       xi(i+1,j+1,k+1,:) = (/real(i),real(j),real(k)/)/geomdim
       if(k==res_z/2) xi(i+1,j+1,k+1,3)= 0.0 ! set highest frequencies to zero
       if(j==res_y/2) xi(i+1,j+1,k+1,2)= 0.0
       if(i==res_x/2) xi(i+1,j+1,k+1,1)= 0.0
 enddo; enddo; enddo

 
 do k = 1, res_z                                       
   do j = 1, res_y
     do i = 1, res_x/2+1
       divergence_field_fft(i,j,k,1) = sum(field_fft(i,j,k,1,:)*xi(i,j,k,:))
       if(vec_tens == 3) then
       divergence_field_fft(i,j,k,2) = sum(field_fft(i,j,k,2,:)*xi(i,j,k,:))
       divergence_field_fft(i,j,k,3) = sum(field_fft(i,j,k,3,:)*xi(i,j,k,:))
       endif
 enddo; enddo; enddo
 divergence_field_fft = divergence_field_fft*img*2.0*pi

 call dfftw_execute_dft_c2r(plan_fft(2), divergence_field_fft, divergence_field)

end subroutine divergence_fft
          
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine divergence(res_x,res_y,res_z,vec_tens,order,geomdim,field,divergence_field)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! calculates divergence field using FDM with variable accuracy
!use vec_tes to decide if tensor (3) or vector (1)
 
 use math       
 implicit none
 integer res_x, res_y, res_z, vec_tens, order
 integer coordinates(6,3)
 real*8 geomdim(3)
 real*8 field(res_x,res_y,res_z,vec_tens,3)
 real*8 divergence_field(res_x,res_y,res_z,vec_tens)
 integer i, j, k, m, l
 real*8, dimension(4,4) :: FDcoefficient = reshape((/ & !from http://en.wikipedia.org/wiki/Finite_difference_coefficients
                                1.0/2.0,      0.0,      0.0,       0.0,&
                                2.0/3.0,-1.0/12.0,      0.0,       0.0,&
                                3.0/4.0,-3.0/20.0,1.0/ 60.0,       0.0,&
                                4.0/5.0,-1.0/ 5.0,4.0/105.0,-1.0/280.0/),&
                               (/4,4/))
 divergence_field = 0.0
 order = order + 1
 do k = 0, res_z-1; do j = 0, res_y-1; do i = 0, res_x-1
   do m = 1, order 
     coordinates(1,:) = mesh_location(mesh_index((/i+m,j,k/),(/res_x,res_y,res_z/)),(/res_x,res_y,res_z/)) + (/1,1,1/) 
     coordinates(2,:) = mesh_location(mesh_index((/i-m,j,k/),(/res_x,res_y,res_z/)),(/res_x,res_y,res_z/)) + (/1,1,1/)
     coordinates(3,:) = mesh_location(mesh_index((/i,j+m,k/),(/res_x,res_y,res_z/)),(/res_x,res_y,res_z/)) + (/1,1,1/) 
     coordinates(4,:) = mesh_location(mesh_index((/i,j-m,k/),(/res_x,res_y,res_z/)),(/res_x,res_y,res_z/)) + (/1,1,1/)
     coordinates(5,:) = mesh_location(mesh_index((/i,j,k+m/),(/res_x,res_y,res_z/)),(/res_x,res_y,res_z/)) + (/1,1,1/) 
     coordinates(6,:) = mesh_location(mesh_index((/i,j,k-m/),(/res_x,res_y,res_z/)),(/res_x,res_y,res_z/)) + (/1,1,1/)
     do l = 1, vec_tens
       divergence_field(i+1,j+1,k+1,l) = divergence_field(i+1,j+1,k+1,l) + FDcoefficient(m,order) * &
                ((field(coordinates(1,1),coordinates(1,2),coordinates(1,3),l,1)- &
                  field(coordinates(2,1),coordinates(2,2),coordinates(2,3),l,1))*real(res_x)/geomdim(1) +&
                 (field(coordinates(3,1),coordinates(3,2),coordinates(3,3),l,2)- &
                  field(coordinates(4,1),coordinates(4,2),coordinates(4,3),l,2))*real(res_y)/geomdim(2) +&
                 (field(coordinates(5,1),coordinates(5,2),coordinates(5,3),l,3)- &
                  field(coordinates(6,1),coordinates(6,2),coordinates(6,3),l,3))*real(res_z)/geomdim(3))  
     enddo
   enddo
 enddo; enddo; enddo
end subroutine divergence


