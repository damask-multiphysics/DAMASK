!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine mesh(res_x,res_y,res_z,geomdim,defgrad_av,centroids,nodes)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 implicit none
 real*8 geomdim(3)
 integer res_x, res_y, res_z
 real*8 wrappedCentroids(res_x+2,res_y+2,res_z+2,3)
 real*8            nodes(res_x+1,res_y+1,res_z+1,3)
 real*8        centroids(res_x  ,res_y  ,res_z  ,3)

 integer, dimension(3,8) :: neighbor = reshape((/ &
                                     0, 0, 0,&
                                     1, 0, 0,&
                                     1, 1, 0,&
                                     0, 1, 0,&
                                     0, 0, 1,&
                                     1, 0, 1,&
                                     1, 1, 1,&
                                     0, 1, 1 &
                                    /), &
                                    (/3,8/))
 
 integer i,j,k,n
 real*8, dimension(3,3) :: defgrad_av
 integer, dimension(3) :: diag, shift, lookup, me, res
 
 diag = 1
 shift = 0
 lookup = 0
 
 res = (/res_x,res_y,res_z/)
 
 wrappedCentroids=0.0
 wrappedCentroids(2:res_x+1,2:res_y+1,2:res_z+1,:) = centroids

 do k=0, res_z+1
   do j=0, res_y+1
    do i=0, res_x+1
       if (k==0 .or. k==res_z+1 .or. &
           j==0 .or. j==res_y+1 .or. &
           i==0 .or. i==res_x+1      ) then
         me = (/i,j,k/)
         shift = sign(abs(res+diag-2*me)/(res+diag),res+diag-2*me)
         lookup = me-diag+shift*res
         wrappedCentroids(i+1,j+1,k+1,:) = centroids(lookup(1)+1,lookup(2)+1,lookup(3)+1,:)- &
                                       matmul(defgrad_av, shift*geomdim)
       endif
 enddo; enddo; enddo
 do k=0, res_z
   do j=0, res_y
    do i=0, res_x
       do n=1,8
         nodes(i+1,j+1,k+1,:) = nodes(i+1,j+1,k+1,:) + wrappedCentroids(i+1+neighbor(n,1),j+1+neighbor(n,2),k+1+neighbor(n,3),:)
 enddo; enddo; enddo; enddo
 nodes = nodes/8.0
  
 end subroutine mesh
 
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine deformed(res_x,res_y,res_z,geomdim,defgrad,defgrad_av,coord_avgCorner)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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