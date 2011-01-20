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
 
!IO.f90 693 2010-11-04 18:18:01Z MPIE\c.kords
!##############################################################
 MODULE IO   
!##############################################################

 CONTAINS
!********************************************************************
! identifies lines without content
!********************************************************************
 pure function IO_isBlank (line)

 use prec, only: pInt
 implicit none

 character(len=*), intent(in) :: line
 character(len=*), parameter :: blank = achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces
 character(len=*), parameter :: comment = achar(35)                               ! comment id '#'
 integer(pInt) posNonBlank, posComment
 logical IO_isBlank
 
 posNonBlank = verify(line,blank)
 posComment  = scan(line,comment)
 IO_isBlank = posNonBlank == 0 .or. posNonBlank == posComment

 return
 
 endfunction

!********************************************************************
! read string value at pos from line
!********************************************************************
 pure function IO_stringValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 character(len=1+positions(pos*2+1)-positions(pos*2)) IO_stringValue

 if (positions(1) < pos) then
   IO_stringValue = ''
 else
   IO_stringValue = line(positions(pos*2):positions(pos*2+1))
 endif
 return

 endfunction
 
!********************************************************************
! read float value at pos from line
!********************************************************************
 pure function IO_floatValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 real(pReal) IO_floatValue

 if (positions(1) < pos) then
   IO_floatValue = 0.0_pReal
 else
   read(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT=*) IO_floatValue
 endif
 return
100 IO_floatValue = huge(1.0_pReal)
 return

 endfunction
 
!********************************************************************
! read int value at pos from line
!********************************************************************
 pure function IO_intValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 integer(pInt) IO_intValue

 if (positions(1) < pos) then
   IO_intValue = 0_pInt
 else
   read(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT=*) IO_intValue
 endif
 return
100 IO_intValue = huge(1_pInt)
 return

 endfunction
 
!********************************************************************
! change character in line to lower case
!********************************************************************
 pure function IO_lc (line)

 use prec, only: pInt
 implicit none

 character (len=*), intent(in) :: line
 character (len=len(line)) IO_lc
 integer(pInt) i

 IO_lc = line
 do i=1,len(line)
    if(64<iachar(line(i:i)) .and. iachar(line(i:i))<91) IO_lc(i:i)=achar(iachar(line(i:i))+32)
 enddo
 return 

 endfunction
 
 !********************************************************************
! locate at most N space-separated parts in line
! return array containing number of parts in line and
! the left/right positions of at most N to be used by IO_xxxVal
!********************************************************************
! pure function IO_stringPos (line,N)
 function IO_stringPos (line,N)

 use prec, only: pReal,pInt
 implicit none

 character(len=*), intent(in) :: line
 character(len=*), parameter :: sep=achar(44)//achar(32)//achar(9)//achar(10)//achar(13) ! comma and whitespaces
 integer(pInt), intent(in) :: N
 integer(pInt)  left,right
 integer(pInt) IO_stringPos(1+N*2)

 IO_stringPos = -1
 IO_stringPos(1) = 0
 right = 0

 do while (verify(line(right+1:),sep)>0)
   left  = right + verify(line(right+1:),sep)
   right = left + scan(line(left:),sep) - 2
   if ( IO_stringPos(1)<N ) then
     IO_stringPos(1+IO_stringPos(1)*2+1) = left
     IO_stringPos(1+IO_stringPos(1)*2+2) = right
   endif
   IO_stringPos(1) = IO_stringPos(1)+1
 enddo

 return

 endfunction
 
END MODULE

program voronoi
 use prec, only: pReal, pInt
 use IO
 implicit none

 logical gotN_Seeds, gotResolution
 logical, dimension(:), allocatable :: grainCheck
 character(len=1024) input_name, output_name, format1, format2, N_Digits, line, key
 integer(pInt) a, b, c, N_Seeds, seedPoint, theGrain, minDistance, myDistance, i, j, k, l, m
 integer(pInt), dimension(3) :: coordinates
 integer(pInt), dimension (1+2*7) ::  posGeom
 real(pReal), dimension(:,:), allocatable :: grainEuler, seeds
 real(pReal), dimension(3) :: dim
 real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal

 print*, '******************************************************************************'
 print*, '                    Spectral Method Problem Set-up'
 print*, '******************************************************************************'
 print*, ''
 print*, 'generates:'
 print*, '    * geom file "_OUTPUT_.geom": Geometrical information for solver'
 print*, '    * material file "material.config": Orientation information for solver'
 print*, '    * "_OUTPUT_.spectral": combined information for solver'
 print*, ''
 write(*, '(A)', advance = 'NO') 'Enter output filename: '
 read(*, *), output_name
 write(*, '(A)', advance = 'NO') 'Enter seed input file (extension .seeds): '
 read(*, *), input_name
  
 open(20, file = trim(input_name)//('.seeds'), status='old', action='read')
 rewind(20)
 do
   read(20,'(a1024)',END = 100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   posGeom = IO_stringPos(line,7)             
   select case ( IO_lc(IO_StringValue(line,posGeom,1)) )
     case ('grains')
       gotN_Seeds = .true.
       N_Seeds = IO_intValue(line,posGeom,2)
     case ('resolution')
       gotResolution = .true.
       do i = 2,6,2
         select case (IO_lc(IO_stringValue(line,posGeom,i)))
           case('a')
             a = IO_intValue(line,posGeom,i+1)
           case('b')
             b = IO_intValue(line,posGeom,i+1)
           case('c')
             c = IO_intValue(line,posGeom,i+1)
         end select
       enddo
   end select
   if (gotN_Seeds .and. gotResolution) exit
 enddo

 100 allocate(grainEuler(N_Seeds,3))
 allocate(Seeds(N_Seeds,3))
 allocate(grainCheck(N_Seeds))
 grainCheck = .false.
 
 print*, 'resolution: ' ,a,b,c
 write(*, '(A)', advance = 'NO') 'New first resolution: '
 read(*, *), a
 write(*, '(A)', advance = 'NO') 'New second resolution: '
 read(*, *), b
 write(*, '(A)', advance = 'NO') 'New third resolution: '
 read(*, *), c
 
 write(*, '(A)', advance = 'NO') 'First dimension: '
 read(*, *), dim(1)
 write(*, '(A)', advance = 'NO') 'Second dimension: '
 read(*, *), dim(2)
 write(*, '(A)', advance = 'NO') 'Third dimension: '
 read(*, *), dim(3)
 
 rewind(20)
 read(20,'(a1024)') line
 posGeom = IO_stringPos(line,2)
 key = IO_stringValue(line,posGeom,2)
 if (IO_lc(key(1:4)) == 'head') then
   do i=1,IO_intValue(line,posGeom,1); read(20,'(a1024)') line; enddo
 else
   rewind(20)
 endif
 do i=1, N_seeds
  read(20,'(a1024)') line
  if (IO_isBlank(line)) cycle                            ! skip empty lines
  posGeom = IO_stringPos(line,6)             
  Seeds(i,1)=IO_floatValue(line,posGeom,1)
  Seeds(i,2)=IO_floatValue(line,posGeom,2)
  Seeds(i,3)=IO_floatValue(line,posGeom,3)
  grainEuler(i,1)=IO_floatValue(line,posGeom,4)
  grainEuler(i,2)=IO_floatValue(line,posGeom,5)
  grainEuler(i,3)=IO_floatValue(line,posGeom,6)
 enddo
 close(20) 
 
 seeds(:,1) = seeds(:,1)*real(a, pReal)
 seeds(:,2) = seeds(:,2)*real(b, pReal)
 seeds(:,3) = seeds(:,3)*real(c, pReal)
 
! calculate No. of digits needed for name of the grains
  i = 1 + int( log10(real( N_Seeds )))
  write(N_Digits, *) i
  N_Digits = adjustl( N_Digits )

!write material.config header and add a microstructure entry for every grain
  open(20, file = trim(output_name)//('_material.config'))
  write(20, '(A)'), '<microstructure>'
  format1 = '(A, I'//trim(N_Digits)//'.'//trim(N_Digits)//', A)'
  format2 = '(A, I'//trim(N_Digits)//', A)'
  do i = 1, N_Seeds
    write(20, trim(format1)), '[Grain', i, ']'
    write(20, '(A)'), 'crystallite 1'
    write(20, trim(format2)), '(constituent)  phase 1   texture ', i, '   fraction 1.0'
  end do

! get random euler angles for every grain, store them in grainEuler and write them to the material.config file
  format2 = '(6(A, F10.6))'
  write(20, '(/, A)'), '<texture>'
    do i = 1, N_Seeds
      write(20, trim(format1)), '[Grain', i, ']'
      write(20, trim(format2)), '(gauss)  phi1 ', grainEuler(i,1), '   Phi ', grainEuler(i,2), &
            &'   Phi2 ', grainEuler(i,3), '   scatter 0   fraction 1'
    end do
  close(20)
  print*, ''
  print*, 'material.config done.'

!write header of geom file
  open(20, file = ((trim(output_name))//'.geom'))
  open(21, file = ((trim(output_name))//'.spectral'))
  write(20, '(A)'), '3 header'
  write(20, '(A, I8, A, I8, A, I8)'), 'resolution  a ', a, '  b ', b, '  c ', c
  write(20, '(A, g15.10, A, g15.10, A, g15.10)'), 'dimension   x ', dim(1), '  y ', dim(2), '  z ', dim(3)
  write(20, '(A)'), 'homogenization  1'

  format1 = '(I'//trim(N_Digits)//'.'//trim(N_Digits)//')'                    ! geom format
  format2 = '(3(tr2, f6.2), 3(I10), I10, a)'                                  ! spectral (Lebensohn) format


! perform voronoi tessellation and write result to files
  do i = 0, a*b*c-1
    minDistance = a*a+b*b+c*c
    do j = 1, N_Seeds
      do k = -1, 1                                                                                            ! left, me, right image
        do l = -1, 1                                                                                          ! front, me, back image
          do m = -1, 1                                                                                        ! lower, me, upper image
            myDistance = ((         mod(i, a)+1 -seeds(j,1)-m*a)**2 + &
                               (mod((i/a), b)+1 -seeds(j,2)-l*b)**2 + &
                           (mod((i/(a*b)), c)+1 -seeds(j,3)-k*c)**2)
            if (myDistance < minDistance) then
              minDistance = myDistance
              theGrain = j
            end if
          end do
        end do
      end do
    end do
    grainCheck(theGrain) = .true.
    write(20, trim(format1)), theGrain
    write(21, trim(format2)), grainEuler(theGrain,1), grainEuler(theGrain,2), grainEuler(theGrain,3), &
                              mod(i, a)+1, mod((i/a), b)+1, mod((i/(a*b)), c)+1, &
                              theGrain, '   1'
  end do
  close(20)
  close(21)
  print*, 'voronoi tesselation done.'
  if (all(grainCheck)) then
    print*, 'all grains mapped!'
  else
    print*, 'only',count(grainCheck),'grains mapped!'
  endif
 
end program voronoi
