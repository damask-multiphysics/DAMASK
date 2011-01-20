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
 integer(pInt) a, b, c, N_Seeds, seedpoint, i
 real(pReal), dimension(:,:), allocatable :: grainEuler, seeds
 real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal
 real(pReal) randomSeed

 print*, '******************************************************************************'
 print*, '                    Voronoi description file'
 print*, '******************************************************************************'
 print*, ''
 print*, 'generates:'
 print*, '    * description file "_GIVEN_NAME_.seeds":'
 print*, ''
 write(*, '(A)', advance = 'NO') 'Please enter value for first resolution: '
 read(*, *), a
 write(*, '(A)', advance = 'NO') 'Please enter value for second resolution: '
 read(*, *), b
 write(*, '(A)', advance = 'NO') 'Please enter value for third resolution: '
 read(*, *), c
 write(*, '(A)', advance = 'NO') 'Please enter No. of Grains: '
 read(*, *), N_Seeds
 write(*, '(A)', advance = 'NO') 'Please enter name of geometry file: '
 read(*, *), filename

 allocate (seedmap(a*b*c)); seedmap = .false.  ! logical to store information which position is occupied by a voronoi seed
 allocate (seeds(N_Seeds,3))
 allocate (grainEuler(N_Seeds,3))
  
 do i=1, N_Seeds
   call random_number(grainEuler(i,1))
   call random_number(grainEuler(i,2))
   call random_number(grainEuler(i,3))
   grainEuler(i,1) = (grainEuler(i,1))*360                     
   grainEuler(i,2) = acos(2.0_pReal*(grainEuler(i,2))-1.0_pReal)*180/pi
   grainEuler(i,3) = grainEuler(i,3)*360
 enddo

!generate random position of seeds for voronoi tessellation
 i = 0
 do while (i /= N_Seeds)
   call random_number(randomSeed)
   seedpoint = int(randomSeed*(a*b*c))
   if (.not.seedmap(seedpoint+1)) then
     seedmap(seedpoint+1) = .true. 
     i = i + 1
     seeds(i,1) = real(mod((seedpoint), a)+1)/real(a, pReal)
     seeds(i,2) = real(mod(((seedpoint)/a), b)+1)/real(b,pReal)
     seeds(i,3) = real(mod(((seedpoint)/(a*b)), c)+1)/real(c,pReal)
   end if
 end do
  
! write description file with orientation and position of each seed
 open(21, file = trim(filename)//('.seeds'))
 write(21, '(A, I2, A, I2, A, I2)'), 'resolution  a ', a, '  b ', b, '  c ', c
 write(21,*), 'grains', N_Seeds
 do i = 1, n_Seeds
   write(21, '(6(F10.6,tr2))'),seeds(i,1), seeds(i,2), seeds(i,3),&
                          grainEuler(i,1), grainEuler(i,2), grainEuler(i,3)              
 end do
 close(21)
end program voronoi
