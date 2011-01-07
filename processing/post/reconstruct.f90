! -*- f90 -*-
subroutine simple(defgrad,res_x,res_y,res_z,geomdimension,current_configuration)
 implicit none
 integer, parameter :: pDouble = selected_real_kind(15,50)
 integer i,j,k
 integer res_x, res_y, res_z
 integer, dimension(3) :: resolution
 real*8, dimension(3) :: geomdimension
 real*8, dimension(3,3) :: temp33_Real
 real*8 current_configuration(res_x, res_y, res_z,3)
 real*8 defgrad(res_x, res_y, res_z,3,3)
 !f2py intent(in) res_x, res_y, res_z
 !f2py intent(in) geomdimension
 !f2py intent(out) current_configuration
 !f2py intent(in) defgrad
 !f2py depend(res_x, res_y, res_z) current_configuration
 !f2py depend(res_x, res_y, res_z) defgrad
 
 resolution(1) = res_x; resolution(2) = res_y; resolution(3) = res_z
 
 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   if((k==1).and.(j==1).and.(i==1)) then
     temp33_Real =0.0_pDouble
   else         
     if((j==1).and.(i==1)) then
       temp33_Real(1,:) = temp33_Real(1,:) + matmul(defgrad(i,j,k,:,:),&
                  (/0.0_pDouble,0.0_pDouble,geomdimension(3)/(real(resolution(3)))/))
       temp33_Real(2,:) = temp33_Real(1,:)
       temp33_Real(3,:) = temp33_Real(1,:)
       current_configuration(i,j,k,:) = temp33_Real(1,:)
     else 
       if(i==1) then
         temp33_Real(2,:) = temp33_Real(2,:) + matmul(defgrad(i,j,k,:,:),&
                   (/0.0_pDouble,geomdimension(2)/(real(resolution(2))),0.0_pDouble/))
         temp33_Real(3,:) = temp33_Real(2,:)
         current_configuration(i,j,k,:) = temp33_Real(2,:)
       else   
         temp33_Real(3,:) = temp33_Real(3,:) + matmul(defgrad(i,j,k,:,:),&
                  (/geomdimension(1)/(real(resolution(1))),0.0_pDouble,0.0_pDouble/))
         current_configuration(i,j,k,:) = temp33_Real(3,:)   
       endif
    endif
  endif
  !print*, current_configuration
 enddo; enddo; enddo
end subroutine simple

subroutine advanced(defgrad,res_x,res_y,res_z,geomdimension,current_configuration)
 implicit none
 integer, parameter :: pDouble = selected_real_kind(15,50)
 integer i,j,k,l
 integer a,b,c
 integer res_x, res_y, res_z
 integer, dimension(3) :: resolution
 real*8, dimension(3) :: geomdimension
 real*8, dimension(3,3) :: temp33_Real, defgrad_av !matrix, stores 3 times a 3dim position vector
 real*8 current_configuration(res_x, res_y, res_z,3)
 real*8 defgrad(res_x, res_y, res_z,3,3)
! some declarations for the wrapping to python 
 !f2py intent(in) res_x, res_y, res_z
 !f2py intent(in) geomdimension
 !f2py intent(out) current_configuration
 !f2py intent(in) defgrad
 !f2py depend(res_x, res_y, res_z) current_configuration
 !f2py depend(res_x, res_y, res_z) defgrad
 do i=1, 3; do j=1,3 !
   defgrad_av(i,j) = sum(defgrad(:,:,:,i,j)) /real(res_x*res_y*res_z)
 enddo; enddo
 
 current_configuration(i,j,k,:) = 0.0_pDouble
 do l=1,6 ! do 6 different paths, by using 3 'master edges'
   select case(l)
     case (1)
       a = 1; b = 2; c = 3
     case (2)
       a = 1; b = 3; c = 2
     case (3)
       a = 2; b = 1; c = 3
     case (4) 
       a = 2; b = 3; c = 1
     case (5)
       a = 3; b = 1; c = 2
     case (6)
       a = 3; b = 2; c = 1
   end select
   
   resolution(a) = res_x; resolution(b) = res_y; resolution(c) = res_z
 
   do k = 1, resolution(c); do j = 1, resolution(b); do i = 1, resolution(a)
     if((k==1).and.(j==1).and.(i==1)) then ! first FP
       temp33_Real = 0.0_pDouble    ! all positions set to zero
     else         
       if((j==1).and.(i==1)) then   ! all FPs on the 'master edge'
         temp33_Real(1,:) = temp33_Real(1,:) + matmul((defgrad(i,j,k-1,:,:)+defgrad(i,j,k-1,:,:))/2.0_pDouble,& !using the average defgrad of the current step
                    (/0.0_pDouble,0.0_pDouble,geomdimension(c)/(real(resolution(c)))/))
         temp33_Real(2,:) = temp33_Real(1,:)
         temp33_Real(3,:) = temp33_Real(1,:)
         current_configuration(i,j,k,:) = current_configuration(i,j,k,:) + temp33_Real(1,:)
       else 
         if(i==1) then
           temp33_Real(2,:) = temp33_Real(2,:) + matmul((defgrad(i,j-1,k,:,:)+defgrad(i,j,k,:,:))/2.0_pDouble,&
                     (/0.0_pDouble,geomdimension(b)/(real(resolution(b))),0.0_pDouble/))
           temp33_Real(3,:) = temp33_Real(2,:)
           current_configuration(i,j,k,:) = current_configuration(i,j,k,:) + temp33_Real(2,:)
         else   
           temp33_Real(3,:) = temp33_Real(3,:) + matmul((defgrad(i-1,j,k,:,:)+defgrad(i,j,k,:,:))/2.0_pDouble,&
                    (/geomdimension(a)/(real(resolution(a))),0.0_pDouble,0.0_pDouble/))
           current_configuration(i,j,k,:) = current_configuration(i,j,k,:) + temp33_Real(3,:)   
         endif
      endif
    endif
 enddo; enddo; enddo ! end of one reconstruction
enddo ! end of 6 reconstructions with different pathes
current_configuration = current_configuration/6.0_pDouble
end subroutine advanced

subroutine mesh(inter,res_x,res_y,res_z,gdim,meshgeom)
 implicit none
 integer, parameter :: pDouble = selected_real_kind(15,50)
 integer i,j,k
 integer res_x, res_y, res_z
 integer, dimension(3) :: resolution
 real*8, dimension(3) :: gdim
 real*8 inter(res_x, res_y, res_z,3)
 real*8 configuration(res_x+2, res_y+2, res_z+2,3)
 real*8 meshgeom(res_x+1, res_y+1, res_z+1,3)
 !f2py intent(in) res_x, res_y, res_z
 !f2py intent(in) inter
 !f2py intent(out) meshgeom
 !f2py intent(in) gdim
 !f2py depend(res_x, res_y, res_z) inter
 !f2py depend(res_x, res_y, res_z) meshgeom
 meshgeom = 0.0_pDouble
 configuration = 0.0_pDouble
 
 resolution(1) = res_x; resolution(2) = res_y; resolution(3) = res_z
 
 
 do k = 1, resolution(3)+1; do j = 1, resolution(2)+1; do i = 1, resolution(1)+1
   if((i==1).and.(j==1).and.(k==1)) then
  configuration(i,j,k,:)=inter(res_x,res_y,res_z,:)-gdim
  configuration(res_x+2,j,k,:)=inter(i,res_y,res_z,:)+(/1.0,-1.0,-1.0/)*gdim
  configuration(i,res_y+2,k,:)=inter(res_x,j,res_z,:)+(/-1.0,1.0,-1.0/)*gdim
  configuration(i,j,res_z+2,:)=inter(res_x,res_y,k,:)+(/-1.0,-1.0,1.0/)*gdim
  configuration(res_x+2,res_y+2,k,:)=inter(i,j,res_z,:)+(/1.0,1.0,-1.0/)*gdim
  configuration(i,res_y+2,res_z+2,:)=inter(res_x,j,k,:)+(/-1.0,1.0,1.0/)*gdim
  configuration(res_x+2,j,res_z+2,:)=inter(i,res_y,k,:)+(/1.0,-1.0,1.0/)*gdim
  configuration(res_x+2,res_y+2,res_z+2,:)=inter(i,j,k,:)+gdim
   endif       
   if((i==1).and.(j==1).and.(k/=1)) then
  configuration(1,1,k,:)=inter(res_x,res_y,k-1,:)+(/-1.0,-1.0,0.0/)*gdim
  configuration(res_x+2,1,k,:)=inter(1,res_y,k-1,:)+(/1.0,-1.0,0.0/)*gdim
  configuration(1,res_y+2,k,:)=inter(res_x,1,k-1,:)+(/-1.0,1.0,0.0/)*gdim 
  configuration(res_x+2,res_y+2,k,:)=inter(1,1,k-1,:)+(/1.0,1.0,0.0/)*gdim
   endif
   if((i==1).and.(j/=1).and.(k==1)) then
  configuration(1,j,1,:)=inter(res_x,j-1,res_z,:)+(/-1.0,0.0,-1.0/)*gdim
  configuration(res_x+2,j,1,:)=inter(1,j-1,res_z,:)+(/1.0,0.0,-1.0/)*gdim
  configuration(1,j,res_z+2,:)=inter(res_x,j-1,1,:)+(/-1.0,0.0,1.0/)*gdim
  configuration(res_x+2,j,res_z+2,:)=inter(1,j-1,1,:)+(/1.0,0.0,1.0/)*gdim
   endif
   if((i/=1).and.(j==1).and.(k==1)) then
     configuration(i,1,1,:)=inter(i-1,res_y,res_z,:)+(/0.0,-1.0,-1.0/)*gdim
     configuration(i,1,res_z+2,:)=inter(i-1,res_y,1,:)+(/0.0,-1.0,1.0/)*gdim
     configuration(i,res_y+2,1,:)=inter(i-1,1,res_z,:)+(/0.0,1.0,-1.0/)*gdim
     configuration(i,res_y+2,res_z+2,:)=inter(i-1,1,1,:)+(/0.0,1.0,1.0/)*gdim
   endif 
   if((i/=1).and.(j/=1).and.(k==1)) then
     configuration(i,j,1,:)=inter(i-1,j-1,res_z,:)+(/0.0,0.0,-1.0/)*gdim
     configuration(i,j,res_z+2,:)=inter(i-1,j-1,1,:)+(/0.0,0.0,1.0/)*gdim
   endif
   if((i==1).and.(j/=1).and.(k/=1)) then
     configuration(1,j,k,:)=inter(res_x,j-1,k-1,:)+(/-1.0,0.0,0.0/)*gdim
     configuration(res_x+2,j,k,:)=inter(i,j-1,k-1,:)+(/1.0,0.0,0.0/)*gdim
   endif
   if((i/=1).and.(j==1).and.(k/=1)) then
     configuration(i,1,k,:)=inter(i-1,res_y,k-1,:)+(/0.0,-1.0,0.0/)*gdim
     configuration(i,res_y+2,k,:)=inter(i-1,1,k-1,:)+(/0.0,1.0,0.0/)*gdim
   endif
   if((i/=1).and.(j/=1).and.(k/=1)) then
     configuration(i,j,k,:)=inter(i-1,j-1,k-1,:)
   endif
 enddo; enddo; enddo

 do k = 1, resolution(3)+1; do j = 1, resolution(2)+1; do i = 1, resolution(1)+1
    meshgeom(i,j,k,:)=((configuration(i,j,k,:)     +configuration(i+1,j,k,:))&
                      +(configuration(i,j+1,k,:)   +configuration(i,j,k+1,:))&
                      +(configuration(i+1,j+1,k,:) +configuration(i,j+1,k+1,:))&
                      +(configuration(i+1,j,k+1,:) +configuration(i+1,j+1,k+1,:)))/8.0_pDouble
 enddo; enddo; enddo
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
