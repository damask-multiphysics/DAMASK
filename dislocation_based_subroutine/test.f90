
!use prec
!use IO
!use constitutive

implicit none

integer, dimension(4) :: array,m

real a,twopi
integer i

array = (/10,2,3,4/)
m = (/1,1,0,0/)

write(*,*) maxval(array,m /= 0)

twopi=2.0*3.1415926
a=1.234

do i=0,-3,-1
  write(*,*) a+i*twopi,modulo(a+i*twopi,twopi)
enddo


!call constitutive_parse_MatTexDat('materials_textures.mpie')

!write(*,*) 'materials_maxN ', materials_maxN
!write(*,*) 'textures_maxN  ', textures_maxN 




end