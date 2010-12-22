! -*- f90 -*-
subroutine simple(defgrad,res_x,res_y,res_z,geomdimension,current_configuration)
 implicit none
!    *** Precision of real and integer variables ***
 integer, parameter :: pReal = selected_real_kind(15)     ! 15 significant digits, up to 1e+-300
 integer, parameter :: pInt  = selected_int_kind(9)           ! up to +- 1e9
 integer(pInt) i,j,k
 integer(pInt) res_x, res_y, res_z
 integer(pInt), dimension(3) :: resolution
 real(pReal) , dimension(3) :: geomdimension
 real(pReal) , dimension(3,3) :: temp33_Real
 real*8 current_configuration(res_x, res_y, res_z,3)
 real*8 defgrad(res_x, res_y, res_z,3,3)
 !f2py intent(in) res_x, res_y, res_z
 !f2py intent(in) geomdimension
 !f2py intent(out) current_configuration
 !f2py intent(in) defgrad
 !f2py depend(res_x, res_y, res_z) current_configuration
 !f2py depend(res_x, res_y, res_z) defgrad
 
 resolution(1) = res_x; resolution(2) = res_y; resolution(3) = res_z
 print*, defgrad
 do k = 1, resolution(3); do j = 1, resolution(2); do i = 1, resolution(1)
   if((k==1).and.(j==1).and.(i==1)) then
     temp33_Real =0.0_pReal
   else         
     if((j==1).and.(i==1)) then
       temp33_Real(1,:) = temp33_Real(1,:) + matmul(defgrad(i,j,k,:,:),&
                  (/0.0_pReal,0.0_pReal,(real(resolution(3), pReal)/geomdimension(3))/))
       temp33_Real(2,:) = temp33_Real(1,:)
       temp33_Real(3,:) = temp33_Real(1,:)
       current_configuration(i,j,k,:) = temp33_Real(1,:)
     else 
       if(i==1) then
         temp33_Real(2,:) = temp33_Real(2,:) + matmul(defgrad(i,j,k,:,:),&
                   (/0.0_pReal,(real(resolution(2),pReal)/geomdimension(2)),0.0_pReal/))
         temp33_Real(3,:) = temp33_Real(2,:)
         current_configuration(i,j,k,:) = temp33_Real(2,:)
       else   
         temp33_Real(3,:) = temp33_Real(3,:) + matmul(defgrad(i,j,k,:,:),&
                  (/(real(resolution(1), pReal)/geomdimension(1)),0.0_pReal,0.0_pReal/))
         current_configuration(i,j,k,:) = temp33_Real(3,:)   
       endif
    endif
  endif  
 enddo; enddo; enddo
end subroutine simple
