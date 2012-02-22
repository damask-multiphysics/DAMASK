!prec.f90 407 2009-08-31 15:09:15Z MPIE\f.roters
!##############################################################
 MODULE prec
!##############################################################

 implicit none
 
!    *** Precision of real and integer variables ***
 integer, parameter :: pReal = selected_real_kind(15,300)     ! 15 significant digits, up to 1e+-300
 integer, parameter :: pInt  = selected_int_kind(9)           ! up to +- 1e9
 integer, parameter :: pLongInt  = 8                          ! should be 64bit

 END MODULE prec

program voronoi
 use prec, only: pReal, pInt
 implicit none

 logical, dimension(:), allocatable :: seedmap
 character(len=1024) filename
 integer(pInt), dimension(3) :: seedcoord
 integer(pInt), dimension(:), allocatable :: rndInit
 integer(pInt) a, b, c, N_Seeds, seedpoint, i, randomSeed, rndSize
 real(pReal), dimension(:,:), allocatable :: grainEuler, seeds
 real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal

 print*, '******************************************************************************'
 print*, '                    Voronoi description file'
 print*, '******************************************************************************'
 print*, '$Id$'
 print*, ''
 print*, 'generates:'
 print*, '    * description file "_OUTPUT_.seeds":'
 print*, ''
 write(*, '(A)', advance = 'NO') 'output seed filename: '
 read(*, *), filename
 write(*, '(A)', advance = 'NO') 'seed of random number generator: '
 read(*, *), randomSeed; randomSeed = max(0_pInt,randomSeed)
 write(*, '(A)', advance = 'NO') 'number of grains: '
 read(*, *), N_Seeds
 write(*, '(A)', advance = 'NO') 'min. Fourier points in x: '
 read(*, *), a
 write(*, '(A)', advance = 'NO') 'min. Fourier points in y: '
 read(*, *), b
 write(*, '(A)', advance = 'NO') 'min. Fourier points in z: '
 read(*, *), c

 allocate (seedmap(a*b*c)); seedmap = .false.  ! logical to store information which position is occupied by a voronoi seed
 allocate (seeds(N_Seeds,3))
 allocate (grainEuler(N_Seeds,3))

 call random_seed(size=rndSize)
 allocate(rndInit(rndSize))
 rndInit = randomSeed
 call random_seed(put=rndInit)
 call random_seed(get=rndInit)
  
 do i=1, N_Seeds
   call random_number(grainEuler(i,1))
   call random_number(grainEuler(i,2))
   call random_number(grainEuler(i,3))
   grainEuler(i,1) = (grainEuler(i,1))*360.0                     
   grainEuler(i,2) = acos(2.0_pReal*(grainEuler(i,2))-1.0_pReal)*180.0/pi
   grainEuler(i,3) = grainEuler(i,3)*360.0
 enddo

!generate random position of seeds for voronoi tessellation
 i = 1
 do while (i <= N_Seeds)
   call random_number(seeds(i,1)); seedcoord(1) = min(a,int(seeds(i,1)*a)+1_pInt)-1_pInt
   call random_number(seeds(i,2)); seedcoord(2) = min(b,int(seeds(i,2)*b)+1_pInt)-1_pInt
   call random_number(seeds(i,3)); seedcoord(3) = min(c,int(seeds(i,3)*c)+1_pInt)-1_pInt
   seedpoint = seedcoord(1) + seedcoord(2)*a + seedcoord(3)*a*b
   if (.not. seedmap(seedpoint+1)) then
     seedmap(seedpoint+1) = .true. 
     i = i + 1
   end if
 end do
  
! write description file with orientation and position of each seed
 open(21, file = trim(filename)//('.seeds'))
 write(21, '(i1,a1,a6)') 4,achar(9),'header'
 write(21, '(A, I8, A, I8, A, I8)') 'resolution  a ', a, '  b ', b, '  c ', c
 write(21, '(A, I8)') 'grains', N_Seeds
 write(21, '(A, I8)') 'random seed ',rndInit(1)
 write(21,'(6(a,a1))') 'x',achar(9),'y',achar(9),'z',achar(9),'phi1',achar(9),'Phi',achar(9),'phi2',achar(9)
 do i = 1, n_Seeds
   write(21, '(6(F10.6,a1))'),seeds(i,1),     achar(9), seeds(i,2),     achar(9), seeds(i,3),     achar(9), &
                              grainEuler(i,1),achar(9), grainEuler(i,2),achar(9), grainEuler(i,3),achar(9)              
 end do
 close(21)
 deallocate (rndInit)
end program voronoi
