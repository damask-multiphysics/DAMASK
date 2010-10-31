program colormap

implicit none
  real startH, endH, startS, endS, startL, endL, h_strich,x,c,m
  integer steps,i,j
  character(len=100) name
  print*, '******************************************************************************'
  print*, '                   write colormap for gmsh'
  print*, '******************************************************************************'
  print*, ''
  write(*, '(A)', advance = 'NO') 'Please enter startvalue for L: '
  read(*, *), startL
  write(*, '(A)', advance = 'NO') 'Please enter endvalue for L: '
  read(*, *), endL
  write(*, '(A)', advance = 'NO') 'Please enter startvalue for S: '
  read(*, *), startS
  write(*, '(A)', advance = 'NO') 'Please enter endvalue for S: '
  read(*, *), endS
  write(*, '(A)', advance = 'NO') 'Please enter steps: '
  read(*, *), steps
  do j=0,360
  write(name, *) j
  name=adjustl(name)
  startH = real(j)
  endH =real(j)
  open(20, file = ('colormap_')//trim(name)//('.map'))
  write(20,*),'View.ColorTable = {'
  if(endH<startH) endH=endH+startH
  do i=0, steps-1
    H_strich = (startH + real(i)*(endH-startH)/real(steps))/60
    if(h_strich>6.0) h_strich = h_strich-6.0
    c = (1- abs(2*(startL + real(i)*(endL-startL)/real(steps))-1))*sqrt(startS + real(i)*(endS-startS)/real(steps))
    x = c*(1- abs(mod(h_strich, real(2))-1))
    m = (startL + real(i)*(endL-startL)/real(steps)) -.5*c
    if ((0.0 <= h_strich).and.(h_strich<1.0)) then
      write(20,*),'{',(c+m)*255,',',(x+m)*255,',',(0.0+m)*255,'},'
    else if ((1.0 <= h_strich).and.(h_strich<2.0)) then
      write(20,*),'{',(x+m)*255,',',(c+m)*255,',',(0.0+m)*255,'},'
    else if ((2.0 <= h_strich).and.(h_strich<3.0)) then
      write(20,*),'{',(0.0+m)*255,',',(c+m)*255,',',(x+m)*255,'},'
    else if ((3.0 <= h_strich).and.(h_strich<4.0)) then
      write(20,*),'{',(0.0+m)*255,',',(x+m)*255,',',(c+m)*255,'},'
    else if ((4.0 <= h_strich).and.(h_strich<5.0)) then
      write(20,*),'{',(x+m)*255,',',(0.0+m)*255,',',(c+m)*255,'},'
    else if ((5.0 <= h_strich).and.(h_strich<=6.0)) then
      write(20,*),'{',(c+m)*255,',',(0.0+m)*255,',',(x+m)*255,'},'  
    endif
  enddo
  H_strich = (startH + real(i)*(endH-startH)/real(steps))/60
  if(h_strich>6.0) h_strich = h_strich-6.0
  c = (1- abs(2*(startL + real(i)*(endL-startL)/real(steps))-1))*(startS + real(i)*(endS-startS)/real(steps))
  x = c*(1- abs(mod(h_strich, real(2))-1))
  m = (startL + real(i)*(endL-startL)/real(steps)) -.5*c
  if ((0.0 <= h_strich).and.(h_strich<1.0)) then
    write(20,*),'{',(c+m)*255,',',(x+m)*255,',',(0.0+m)*255,'}'
  else if ((1.0 <= h_strich).and.(h_strich<2.0)) then
    write(20,*),'{',(x+m)*255,',',(c+m)*255,',',(0.0+m)*255,'}'
  else if ((2.0 <= h_strich).and.(h_strich<3.0)) then
    write(20,*),'{',(0.0+m)*255,',',(c+m)*255,',',(x+m)*255,'}'
  else if ((3.0 <= h_strich).and.(h_strich<4.0)) then
    write(20,*),'{',(0.0+m)*255,',',(x+m)*255,',',c+m,'}'
  else if ((4.0 <= h_strich).and.(h_strich<5.0)) then
    write(20,*),'{',(x+m)*255,',',(0.0+m)*255,',',(c+m)*255,'}'
  else if ((5.0 <= h_strich).and.(h_strich<=6.0)) then
    write(20,*),'{',(c+m)*255,',',(0.0+m)*255,',',(x+m)*255,'}'  
  endif
  write(20,*),'};'
  enddo
end program 