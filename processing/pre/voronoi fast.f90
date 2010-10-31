program voronoi
  use prec, only: pReal, pInt

  implicit none

  logical, dimension(:), allocatable :: seedmap
  character(len=1024) name, format, format2, N_Digits
  character choice
  integer(pInt) a, b, c, ab, abc, abc_Red, N_Seeds, seedPoint, minDistance, myDistance, i, j, k, l, m
  integer(pInt), dimension(:), allocatable :: seeds, grainMap, visual_Case
  integer(pInt) coordinates(3)
  real(pReal), dimension(:), allocatable :: grainEuler
  real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal
  real(pReal) randomSeed

  print*, '******************************************************************************'
  print*, '                    Spectral Method Problem Set-up'
  print*, '******************************************************************************'
  print*, ''
  print*, 'generates:'
  print*, '    * mesh file "_GIVEN_NAME_.mesh": Geometrical information for solver'
  print*, '    * material file "material.config": Orientation information for solver'
  print*, '    * "_GIVEN_NAME_.spectral": combined information for solver'
  print*, 'optional output:'
  print*, '    * view file "_GIVEN_NAME_2D.msh": Information for visualization in gmsh'
  print*, '    * view file "_GIVEN_NAME_3D.msh": Information for visualization in gmsh'
  print*, ''
  print*, 'hints:'
  print*, '    * a+b+c should not exeed 30'
  print*, '    * file extension is added to given name'
  print*, ''
  write(*, '(A)', advance = 'NO') 'Please enter value for a: '
  read(*, *), a
  write(*, '(A)', advance = 'NO') 'Please enter value for b: '
  read(*, *), b
  write(*, '(A)', advance = 'NO') 'Please enter value for c: '
  read(*, *), c
  write(*, '(A)', advance = 'NO') 'Please enter No. of Grains: '
  read(*, *), N_Seeds
  write(*, '(A)', advance = 'NO') 'Please enter name of mesh file: '
  read(*, *), name
  write(*, '(A)', advance = 'NO') 'Should the visualization files be generated (y/n)? '
  read(*, *), choice

! calculate No. of digits needed for name of the grains
  i = 1 + int( log10(real( N_Seeds )))
  write(N_Digits, *) i
  N_Digits = adjustl( N_Digits )
  allocate(grainEuler(3*N_Seeds))
  
!write material.config header and add a microstructure entry for every grain
  open(20, file = trim(name)//('_material.config'))
  write(20, '(A)'), '<microstructure>'
  format = '(A, I'//trim(N_Digits)//'.'//trim(N_Digits)//', A)'
  format2 = '(A, I'//trim(N_Digits)//', A)'
  do i = 1, N_Seeds
    write(20, trim(format)), '[Grain', i, ']'
    write(20, '(A)'), 'crystallite	1'
    write(20, trim(format2)), '(constituent)  phase 1   texture ', i, '   fraction 1.0'
  end do

! get random euler angles for every grain, store them in grainEuler and write them to the material.config file
  format2 = '(A, F10.6, A, F10.6, A, F10.6, A)'
  write(20, '(/, A)'), '<texture>'
    do i = 1, N_Seeds
      call random_number(grainEuler(i*3 -2))
      call random_number(grainEuler(i*3 -1))
      call random_number(grainEuler(i*3 -0))
      grainEuler(i*3 -2) = (grainEuler(i*3 -2))*360
      grainEuler(i*3 -1) = acos(2.0_pReal*(grainEuler(i*3 -1))-1.0_pReal)*180/pi
      grainEuler(i*3 -0) = grainEuler(i*3)*360
      write(20, trim(format)), '[Grain', i, ']'
      write(20, trim(format2)), '(gauss)  phi1 ', grainEuler(i*3-2), '   Phi ', grainEuler(i*3-1), &
            &'   Phi2 ', grainEuler(i*3), '   scatter 0   fraction 1'
    end do
  close(20)
  print*, ''
  print*, 'material config file is written out'

!write header of mesh file, should be done before the following change of variables
  open(20, file = ((trim(name))//'.mesh'))
  write(20, '(A, I2, A, I2, A, I2)'), 'resolution  a ', 2**a, '  b ', 2**b, '  c ', 2**c
  write(20, '(A, I4, A, I4, A, I4)'), 'dimension   x ', 2**a, '  y ', 2**b, '  z ', 2**c
  write(20, '(A)'), 'homogenization  1'

!initialize varibles, change values of some numbers for faster execution
  a = 2**a
  b = 2**b
  c = 2**c
  ab = a * b
  abc = a * b * c
  abc_Red = abc -(a-1)*(b-1)*(c-1)
  format = '(I'//trim(N_Digits)//'.'//trim(N_Digits)//')'

  allocate (seedmap(abc)); seedmap = .false.
  allocate (seeds(3*N_Seeds))
  allocate (grainMap(abc))
  allocate (visual_Case(abc_Red))
  k = 1

!build array with x-y-z-coordinates of each point
  do i = 1, abc
    coordinates(1) = mod((i-1), a) +1
    coordinates(2) = mod(((i-1)/a), b) +1
    coordinates(3) = mod(((i-1)/(ab)), c) +1
    if((coordinates(3) == 1)) then
      visual_Case(k) = i
      k = k +1
    else
      if((coordinates(2) == 1)) then
        visual_Case(k) = i
        k = k +1
      else
        if((coordinates(1) == 1)) then
          visual_Case(k) = i
          k = k +1
        else
        end if
      end if
    end if
  end do

!generate random position of seeds for voronoi tessellation
  i = 0
  do while (i /= N_Seeds)
    call random_number(randomSeed)
    seedpoint = int(randomSeed*(abc))
    if (.not.seedmap(seedpoint+1)) then
      seedmap(seedpoint+1) = .true.
      seeds(i*3+1) = (mod((seedpoint), a)+1)
      seeds(i*3+2) = (mod(((seedpoint)/a), b)+1)
      seeds(i*3+3) = (mod(((seedpoint)/(ab)), c)+1)
      i = i +1
    else
    end if
  end do

! perform voronoi tessellation and write result to file and to grainMap
  do i = 1, abc
    minDistance = a*a+b*b+c*c
    do j = 1, N_Seeds
      do k = -1, 1
        do l = -1, 1
          do m = -1, 1
            myDistance = ((mod((i-1), a) +1-seeds(j*3-2)+m*a)**2+&
                     (mod(((i-1)/a), b) +1-seeds(j*3-1)+l*b)**2+&
                     (mod(((i-1)/(ab)), c) +1-seeds(j*3-0)+k*c)**2)
            if (myDistance < minDistance) then
              minDistance = myDistance
              grainMap(i) = j
            end if
          end do
        end do
      end do
    end do
    write(20, format), grainMap(i)
  end do
  close(20)
  print*, 'voronoi tesselation finished'

  open(20, file = ((trim(name))//'.spectral'))
  format = '(tr2, f6.2, tr2, f6.2, tr2, f6.2, I10, I10, I10, I10, a)'
  do i = 1, abc
    j = grainMap(i)
    write(20, trim(format)), grainEuler(j*3-2), grainEuler(j*3-1), grainEuler(j*3), &
	                       &mod((i-1), a)+1, mod(((i-1)/a), b)+1, mod(((i-1)/(ab)), c)+1, &
	                       &j, '   1'
  end do
  print*, 'mesh files are written out'
  
!write visualization files (in case wanted)
  if (choice == 'y' .or. choice == 'Y') then
    print*, 'for more information on gmsh: http://geuz.org/gmsh/'
	
! write full mesh out
    open(20, file = ((trim(name))//'_3D.msh'))
    write(20, '(A, /, A, /, A, /, A, /, I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', abc_Red
    do j = 1, abc_Red
      i = visual_Case(j)
      write(20, '(I10, I10, I10, I10)'), i, mod((i-1), a) +1, mod(((i-1)/a), b) +1, mod(((i-1)/(ab)), c) +1
    end do
    write(20, '(A, /, A, /, I10)'), '$EndNodes', '$Elements', abc_Red
    do j = 1, abc_Red
      i = visual_case(j)
      write(20, '(I10, A, I10, A, I10)'), i, ' 15 2', grainMap(i), ' 2', i
    end do
    write(20, '(A)'), '$EndElements'
    write(20, '(A, /, A, /, A, /, A, /, A, /, A, /, A, /, A, /, I10)'), '$NodeData', '1', '"Grain No."', '1', &
                                                                 &'0.0', '3', '0', '1', abc_Red
    do j = 1, abc_Red
      i = visual_case(j)
      write(20, '(I10, tr2, I10)'), i, grainMap(i)
    end do
    write(20, *), '$EndNodeData'
    close(20)

! write 2d mesh out
    open(20, file = ((trim(name))//'_2D.msh'))
    write(20, '(A, /, A, /, A, /, A, /, I10)'), '$MeshFormat', '2.1 0 8', '$EndMeshFormat', '$Nodes', ab
    do i = 1, ab
      write(20, '(I10, I10, I10, I10)'), i, mod((i-1), a) +1, mod(((i-1)/a), b) +1, mod(((i-1)/(ab)), c) +1
    end do
    write(20, '(A, /, A, /, I10)'), '$EndNodes', '$Elements', ab
    do j = 1, ab
      write(20, '(I10, A, I10, A, I10)'), j, ' 15 2', grainMap(j), ' 2', j
    end do
    write(20, '(A)'), '$EndElements'
    write(20, '(A, /, A, /, A, /, A, /, A, /, A, /, A, /, A, /, I10)'), '$NodeData', '1', '"Grain No."', &
                                                    &'1', '0.0', '3', '0', '1', ab
    do j = 1, ab
      write(20, '(I10, tr2, I10)'), j, grainMap(j)
    end do
    write(20, *), '$EndNodeData'
    close(20)
    print*, 'visualization files are written out'
    else
  end if
end program voronoi
