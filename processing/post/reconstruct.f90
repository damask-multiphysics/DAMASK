! -*- f90 -*-
function coordinates2(res_x,res_y,res_z,geomdimension,defgrad)
 implicit none
 integer, parameter :: pDouble = selected_real_kind(15,50)
 integer i,j,k, l,m, s,o, loop, res_x, res_y, res_z
 integer, dimension(3) :: res, init, oppo, me, rear
 real*8, dimension(3,3) :: defgrad_av
 integer, dimension(3) :: resolution
 real*8, dimension(3) :: geomdimension, myStep
 real*8, dimension(3,3) :: temp33_Real
 real*8, dimension(3,res_x, res_y, res_z) :: coordinates2
 real*8, dimension(3,3,res_x, res_y, res_z) :: defgrad
 real*8, dimension(3,res_x,res_y,res_z,8) :: cornerCoords
 real*8, dimension(3,2+res_x,2+res_y,2+res_z,6,8) :: coord
 !f2py intent(in) res_x
 !f2py intent(in) res_y
 !f2py intent(in) res_z
 !f2py intent(in) geomdimension
 !f2py intent(in) defgrad
 !f2py intent(out) coordinates2
 !f2py depend(res_x, res_y, res_z) coordinates2
 !f2py depend(res_x, res_y, res_z) defgrad
 
 ! integer, dimension(3,8) :: corner = reshape((/ &
                                     ! 0, 0, 0,&
                                     ! 1, 0, 0,&
                                     ! 1, 1, 0,&
                                     ! 0, 1, 0,&
                                     ! 1, 1, 1,&
                                     ! 0, 1, 1,&
                                     ! 0, 0, 1,&
                                     ! 1, 0, 1 &
                                    ! /), &
                                    ! (/3,8/))
 
 ! integer, dimension(3,8) :: step = reshape((/ &
                                     ! 1, 1, 1,&
                                    ! -1, 1, 1,&
                                    ! -1,-1, 1,&
                                     ! 1,-1, 1,&
                                    ! -1,-1,-1,&
                                     ! 1,-1,-1,&
                                     ! 1, 1,-1,&
                                    ! -1, 1,-1 &
                                    ! /), &
                                    ! (/3,8/))
 
 ! integer, dimension(3,6) :: order = reshape((/ &
                                     ! 1, 2, 3,&
                                     ! 1, 3, 2,&
                                     ! 2, 1, 3,&
                                     ! 2, 3, 1,&
                                     ! 3, 1, 2,&
                                     ! 3, 2, 1 &
                                    ! /), &
                                    ! (/3,6/))
 
 resolution = (/res_x,res_y,res_z/)
 
 write(6,*) 'defgrad', defgrad
 
 do i=1, 3; do j=1,3
   defgrad_av(i,j) = sum(defgrad(i,j,:,:,:)) /real(res_x*res_y*res_z)
 enddo; enddo
 
 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   if((k==1).and.(j==1).and.(i==1)) then
     temp33_Real = real(0.0)
   else         
     if((j==1).and.(i==1)) then
       temp33_Real(1,:) = temp33_Real(1,:) + matmul(defgrad(:,:,i,j,k),&
                  (/real(0.0),real(0.0),real(geomdimension(3))/(real(resolution(3)))/))
       temp33_Real(2,:) = temp33_Real(1,:)
       temp33_Real(3,:) = temp33_Real(1,:)
       coordinates2(:,i,j,k) = temp33_Real(1,:)
     else 
       if(i==1) then
         temp33_Real(2,:) = temp33_Real(2,:) + matmul(defgrad(:,:,i,j,k),&
                   (/real(0.0),real(geomdimension(2))/(real(resolution(2))),real(0.0)/))
         temp33_Real(3,:) = temp33_Real(2,:)
         coordinates2(:,i,j,k) = temp33_Real(2,:)
       else   
         temp33_Real(3,:) = temp33_Real(3,:) + matmul(defgrad(:,:,i,j,k),&
                  (/real(geomdimension(1))/(real(resolution(1))),real(0.0),real(0.0)/))
         coordinates2(:,i,j,k) = temp33_Real(3,:)   
       endif
    endif
  endif
 enddo; enddo; enddo
  do i=1, res_x; do j = 1, res_y; do k = 1, res_z  
   coordinates2(:,i,j,k) = coordinates2(:,i,j,k)+ matmul(defgrad_av,(/geomdimension(1)/real(res_x),geomdimension(2)/real(res_y),geomdimension(3)/real(res_z)/))
 enddo; enddo; enddo
 
 res = (/res_x,res_y,res_z/)
 do i=1,3; do j=1,3
   defgrad_av(i,j) = sum(defgrad(i,j,:,:,:)) /real(res(1)*res(2)*res(3))
 enddo; enddo

 ! do s = 1,8
   ! init = corner(:,s)*(res-(/1,1,1/))
   ! oppo = corner(:,1+mod(s-1+4,8))*(res-(/1,1,1/))
     ! do o = 1,6
       ! do k = init(order(3,o)),oppo(order(3,o)),step(order(3,o),s)
         ! rear(order(2,o)) = init(order(2,o))
          ! do j = init(order(2,o)),oppo(order(2,o)),step(order(2,o),s)
          ! rear(order(1,o)) = init(order(1,o))
           ! do i = init(order(1,o)),oppo(order(1,o)),step(order(1,o),s)
 ! !            print*, order(:,o)
             ! me(order(1,o)) = i
             ! me(order(2,o)) = j
             ! me(order(3,o)) = k
               ! ! print*, me
! !             if (all(me == init)) then
! !               coord(:,1+me(1),1+me(2),1+me(3),o,s) = 0.0 !&
                 ! ! geomdimension*(matmul(defgrad_av,real(corner(:,s),pDouble)) + &
                                ! ! matmul(defGrad(:,:,1+me(1),1+me(2),1+me(3)),step(:,s)/res/2.0_pDouble))
             ! ! else
             ! ! myStep = (me-rear)*geomdimension/res
             ! ! coord(:,1+me(1),1+me(2),1+me(3),o,s) = coord(:,1+rear(1),1+rear(2),1+rear(3),o,s) + &
                                                    ! ! 0.5_pDouble*matmul(defGrad(:,:,1+me(1),1+me(2),1+me(3)) + &
                                                                       ! ! defGrad(:,:,1+rear(1),1+rear(2),1+rear(3)),&
                                                                       ! ! myStep)
             ! ! endif
           ! ! rear = me
           ! enddo
          ! enddo         
       ! enddo         
     ! enddo                                             ! orders
   ! ! cornerCoords(:,:,:,:,s) = sum(coord(:,:,:,:,:,s),5)/6.0_pDouble
 ! enddo                                               ! corners
 
 ! coordinates = sum(cornerCoords(:,:,:,:,:),5)/8.0_pDouble ! plain average no shape functions...
end function


 
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
subroutine coordinates5(res_x,res_y,res_z,geomdimension,defgrad)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 implicit none
 integer, parameter :: pDouble = selected_real_kind(15,50)
 integer i,j,k, l,m, s,o, loop, res_x, res_y, res_z
 integer, dimension(3) :: res, init, oppo, me, rear
 real*8, dimension(3) :: geomdimension
 real*8, dimension(3) :: myStep
 real*8, dimension(3,3) :: defGrad_av
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
                                    
 real*8 defGrad(3,3,res_x,res_y,res_z)
 real*8 coordinates(3,res_x,res_y,res_z) 
 real*8, dimension(3,res_x,res_y,res_z,8) :: cornerCoords
 real*8, dimension(3,2+res_x,2+res_y,2+res_z,6,8) :: coord
 !f2py intent(in) res_x, res_y, res_z
 !f2py intent(in) geomdimension
 !f2py intent(out) coordinates
 !f2py intent(in) defgrad
 !f2py depend(res_x, res_y, res_z) coordinates
 !f2py depend(res_x, res_y, res_z) defgrad 
 !f2py depend(res_x, res_y, res_z) cornerCoords
 !f2py depend(res_x, res_y, res_z) coord

 res = (/res_x,res_y,res_z/)
 do i=1,3; do j=1,3
   defgrad_av(i,j) = sum(defgrad(i,j,:,:,:)) /real(res(1)*res(2)*res(3))
 enddo; enddo
 print*, 'defgra', defgrad
 do s = 1,8
   init = corner(:,s)*(res-(/1,1,1/))
   oppo = corner(:,1+mod(s-1+4,8))*(res-(/1,1,1/))
   do o = 1,6
     do k = init(order(3,o)),oppo(order(3,o)),step(order(3,o),s)
       rear(order(2,o)) = init(order(2,o))
       do j = init(order(2,o)),oppo(order(2,o)),step(order(2,o),s)
         rear(order(1,o)) = init(order(1,o))
         do i = init(order(1,o)),oppo(order(1,o)),step(order(1,o),s)
           me(order(1,o)) = i
           me(order(2,o)) = j
           me(order(3,o)) = k
           if (all(me == init)) then
             coord(:,1+me(1),1+me(2),1+me(3),o,s) = &
               geomdimension*(matmul(defgrad_av,real(corner(:,s),pDouble)) + &
                              matmul(defGrad(:,:,1+me(1),1+me(2),1+me(3)),step(:,s)/res/2.0_pDouble))
           else
             myStep = (me-rear)*geomdimension/res
             coord(:,1+me(1),1+me(2),1+me(3),o,s) = coord(:,1+rear(1),1+rear(2),1+rear(3),o,s) + &
                                                    0.5_pDouble*matmul(defGrad(:,:,1+me(1),1+me(2),1+me(3)) + &
                                                                       defGrad(:,:,1+rear(1),1+rear(2),1+rear(3)),&
                                                                       myStep)
           endif
           rear = me
         enddo
       enddo         
     enddo         
   enddo                                             ! orders
   cornerCoords(:,:,:,:,s) = sum(coord(:,:,:,:,:,s),5)/6.0_pDouble
 enddo                                               ! corners
 
 coordinates = sum(cornerCoords(:,:,:,:,:),5)/8.0_pDouble ! plain average no shape functions...

end subroutine

! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
subroutine mesh(res_x,res_y,res_z,geomdim,defgrad_av,centroids,nodes)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
 implicit none
 integer, parameter :: pDouble = selected_real_kind(15,50)
 integer i,j,k, n
 integer res_x, res_y, res_z
 integer, dimension(3) :: res, shift, lookup, me
 real*8, dimension(3) :: geomdim, diag = (/1,1,1/)
 real*8, dimension(3,3) :: defgrad_av
 real*8, dimension(3,res_x, res_y, res_z) :: centroids
 real*8, dimension(3,res_x+2, res_y+2, res_z+2) :: wrappedCentroids
 real*8, dimension(3,res_x+1, res_y+1, res_z+1) :: nodes
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
 !f2py intent(in) res_x, res_y, res_z
 !f2py intent(in) centroids
 !f2py intent(in) defgrad_av
 !f2py intent(in) geomdim
 !f2py intent(out) nodes
 !f2py depend(res_x, res_y, res_z) centroids
 !f2py depend(res_x, res_y, res_z) nodes

 res = (/res_x,res_y,res_z/)
 wrappedCentroids(:,2:res(1)+1,2:res(2)+1,2:res(3)+1) = centroids
 
 do k = 0,res(3)+1
   do j = 0,res(2)+1
     do i = 0,res(1)+1
       if (&
         k==0 .or. k==res(3)+1 .or. &
         j==0 .or. j==res(2)+1 .or. &
         i==0 .or. i==res(1)+1      &
         ) then
         me = (/i,j,k/)
         shift = (res+diag-2*me)/(res+diag)
         lookup = me+shift*res
         wrappedCentroids(:,1+i,1+j,1+k) = centroids(:,lookup(1),lookup(2),lookup(3)) - &
                                           matmul(defgrad_av,shift * geomdim)
       endif
     enddo
   enddo
 enddo
 
 nodes = 0.0_pDouble
 
 do k = 0,res(3)
   do j = 0,res(2)
     do i = 0,res(1)
       
       do n = 1,8
         nodes(:,1+i,1+j,1+k) = nodes(:,1+i,1+j,1+k) + wrappedCentroids(:,1+i+neighbor(1,n),&
                                                                          1+j+neighbor(2,n),&
                                                                          1+k+neighbor(3,n))
       enddo
       nodes(:,1+i,1+j,1+k) = nodes(:,1+i,1+j,1+k) / 8.0_pDouble

     enddo 
   enddo 
 enddo 

end subroutine mesh


!below some code I used for gmsh postprocessing. Might be helpful
!!gmsh output
 ! character(len=1024) :: nriter
 ! character(len=1024) :: nrstep
 ! character(len=1024) :: nrloadcase
 ! real(pReal), dimension(:,:,:,:), allocatable ::    displacement
 ! real(pReal), dimension(3,3) ::                     temp33_Real2
 ! real(pReal), dimension(3,3,3) ::                   Eigenvectorbasis
 ! real(pReal), dimension(3) ::                       Eigenvalue
 ! real(pReal) determinant
!!gmsh output
!!Postprocessing (gsmh output)
     ! do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       ! if((k==1).and.(j==1).and.(i==1)) then
           ! temp33_Real =0.0_pReal
       ! else         
         ! if((j==1).and.(i==1)) then
           ! temp33_Real(1,:) = temp33_Real(1,:) + math_mul33x3(defgrad(i,j,k,:,:),&
                            ! (/0.0_pReal,0.0_pReal,(real(resolution(3))/meshdimension(3))/))
           ! temp33_Real(2,:) = temp33_Real(1,:)
           ! temp33_Real(3,:) = temp33_Real(1,:)
           ! displacement(i,j,k,:) = temp33_Real(1,:)
         ! else 
           ! if(i==1) then
             ! temp33_Real(2,:) = temp33_Real(2,:) + math_mul33x3(defgrad(i,j,k,:,:),&
                       ! (/0.0_pReal,(real(resolution(2))/meshdimension(2)),0.0_pReal/))
             ! temp33_Real(3,:) = temp33_Real(2,:)
             ! displacement(i,j,k,:) = temp33_Real(2,:)
           ! else   
             ! temp33_Real(3,:) = temp33_Real(3,:) + math_mul33x3(defgrad(i,j,k,:,:),&
                       ! (/(real(resolution(1))/meshdimension(1)),0.0_pReal,0.0_pReal/))
             ! displacement(i,j,k,:) = temp33_Real(3,:)   
           ! endif
         ! endif
       ! endif  
     ! enddo; enddo; enddo

     ! write(nrloadcase, *) loadcase; write(nriter, *) iter; write(nrstep, *) steps
     ! open(589,file = 'stress' //trim(adjustl(nrloadcase))//'-'//trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh')
     ! open(588,file = 'logstrain'//trim(adjustl(nrloadcase))//'-'//trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh')
     ! write(589, '(4(A, /), I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn
     ! write(588, '(4(A, /), I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', prodnn

     ! ielem = 0_pInt
     ! do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
       ! ielem = ielem + 1
       ! write(589, '(I10, 3(tr2, E12.6))'), ielem, displacement(i,j,k,:) !for deformed configuration
       ! write(588, '(I10, 3(tr2, E12.6))'), ielem, displacement(i,j,k,:)        
       !! write(589, '(4(I10,tr2))'), ielem, i-1,j-1,k-1 !for undeformed configuration
       !!write(588, '(4(I10,tr2))'), ielem, i-1,j-1,k-1 
     ! enddo; enddo; enddo

     ! write(589, '(2(A, /), I10)'), '$EndNodes', '$Elements', prodnn
     ! write(588, '(2(A, /), I10)'), '$EndNodes', '$Elements', prodnn

     ! do i = 1, prodnn
       ! write(589, '(I10, A, I10)'), i, ' 15 2    1     2', i
       ! write(588, '(I10, A, I10)'), i, ' 15 2    1     2', i
     ! enddo

     ! write(589, '(A)'), '$EndElements'
     ! write(588, '(A)'), '$EndElements'
     ! write(589, '(8(A, /), I10)'), '$NodeData', '1','"'//trim(adjustl('stress'//trim(adjustl(nrloadcase))//'-'//&
                            ! trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh'))//'"','1','0.0', '3', '0', '9', prodnn
     ! write(588, '(8(A, /), I10)'), '$NodeData', '1','"'//trim(adjustl('logstrain'//trim(adjustl(nrloadcase))//'-'//&
                            ! trim(adjustl(nrstep))//'-'//trim(adjustl(nriter))//'_cpfem.msh'))//'"','1','0.0', '3', '0', '9', prodnn
     ! ielem = 0_pInt
     ! do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
            ! ielem = ielem + 1
             ! write(589, '(i10, 9(tr2, E14.8))'), ielem, cstress_field(i,j,k,:,:)
             ! call math_pDecomposition(defgrad(i,j,k,:,:),temp33_Real2,temp33_Real,errmatinv)  !store R in temp33_Real
              ! call math_invert3x3(temp33_Real, temp33_Real2, determinant, errmatinv)           !inverse of R in temp33_Real2
             ! temp33_Real = math_mul33x33(defgrad(i,j,k,:,:), temp33_Real2)                    ! v = F o inv(R), store in temp33_Real
             ! call math_spectral1(temp33_Real,Eigenvalue(1), Eigenvalue(2), Eigenvalue(3),&
                              ! Eigenvectorbasis(1,:,:),Eigenvectorbasis(2,:,:),Eigenvectorbasis(3,:,:))
             ! eigenvalue = log(sqrt(eigenvalue))
             ! temp33_Real = eigenvalue(1)*Eigenvectorbasis(1,:,:)+eigenvalue(2)*Eigenvectorbasis(2,:,:)+eigenvalue(3)*Eigenvectorbasis(3,:,:)
             ! write(588, '(i10, 9(tr2, E14.8))'), ielem, temp33_Real
     ! enddo; enddo; enddo

     ! write(589, *), '$EndNodeData'
     ! write(588, *), '$EndNodeData'
     ! close(589); close(588); close(540) 
   ! enddo  ! end looping over steps in current loadcase
 ! enddo    ! end looping over loadcases
! close(539); close(538)



! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
function coordinates3(res_x,res_y,res_z,geomdimension,defgrad)
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 implicit none
 integer i,j,k, l,m, s,o, loop, res_x, res_y, res_z
 real*8 defgrad_av(3,3)
 real*8 geomdimension(3)
 integer res(3)
 real*8 defgrad(3,3,res_x,res_y,res_z)
 real*8 coordinates3(3,res_x,res_y,res_z) 

 !f2py intent(in) res_x, res_y, res_z
 !f2py depend(res_x, res_y, res_z) coordinates3
 !f2py depend(res_x, res_y, res_z) defgrad 

 !f2py intent(in) geomdimension
 !f2py intent(out) coordinates3
 !f2py intent(in) defgrad
 
 res = (/res_x,res_y,res_z/)
 do i=1,3; do j=1,3
   defgrad_av(i,j) = sum(defgrad(i,j,:,:,:)) /real(res(1)*res(2)*res(3))
 enddo; enddo
 print*, 'defgra', defgrad
 coordinates3 = 0.0
end function