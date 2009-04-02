
!##############################################################
 MODULE IO   
!##############################################################

 CONTAINS
!---------------------------
! function IO_open_file(unit,relPath)
! function IO_open_inputFile(unit)
! function IO_hybridIA(Nast,ODFfileName)
! private function hybridIA_reps(dV_V,steps,C)
! function IO_stringPos(line,maxN)
! function IO_stringValue(line,positions,pos)
! function IO_floatValue(line,positions,pos)
! function IO_intValue(line,positions,pos)
! function IO_fixedStringValue(line,ends,pos)
! function IO_fixedFloatValue(line,ends,pos)
! function IO_fixedFloatNoEValue(line,ends,pos)
! function IO_fixedIntValue(line,ends,pos)
! function IO_continousIntValues(unit,maxN)
! function IO_lc(line)
! subroutine IO_lcInplace(line)
! subroutine IO_error(ID)
! subroutine IO_warning(ID)
!---------------------------



!********************************************************************
! open existing file to given unit
! path to file is relative to working directory
!********************************************************************
 logical FUNCTION IO_open_file(unit,relPath)

 use prec, only: pInt
 implicit none

 character(len=*), parameter :: pathSep = achar(47)//achar(92) ! /, \
 character(len=*) relPath
 integer(pInt) unit
 character(256) path

 inquire(6, name=path) ! determine outputfile
 open(unit,status='old',err=100,file=path(1:scan(path,pathSep,back=.true.))//relPath)
 IO_open_file = .true.
 return
100 IO_open_file = .false.
 return
 
 END FUNCTION


!********************************************************************
! open FEM inputfile to given unit
!********************************************************************
 logical FUNCTION IO_open_inputFile(unit)

 use prec, only: pReal, pInt
 implicit none

 character(256) outName
 integer(pInt) unit, extPos
 character(3) ext

 inquire(6, name=outName) ! determine outputfileName
 extPos = len_trim(outName)-2
 if(outName(extPos:extPos+2)=='out') then
     ext='dat' ! MARC
 else
     ext='inp' ! ABAQUS
 end if
 open(unit,status='old',err=100,file=outName(1:extPos-1)//ext)
 IO_open_inputFile = .true.
 return
100 IO_open_inputFile = .false.
 return

 END FUNCTION


!********************************************************************
! hybrid IA repetition counter
!********************************************************************
 FUNCTION hybridIA_reps(dV_V,steps,C)

 use prec, only: pReal, pInt
 implicit none
 
 integer(pInt), intent(in), dimension(3) :: steps
 integer(pInt) hybridIA_reps, phi1,Phi,phi2
 real(pReal), intent(in), dimension(steps(3),steps(2),steps(1)) :: dV_V
 real(pReal), intent(in) :: C
 
 hybridIA_reps = 0_pInt
 do phi1=1,steps(1)
   do Phi =1,steps(2)
     do phi2=1,steps(3)
       hybridIA_reps = hybridIA_reps+nint(C*dV_V(phi2,Phi,phi1), pInt)
     end do
   end do
 end do
 return
 
 END FUNCTION


!********************************************************************
! hybrid IA sampling of ODFfile
!********************************************************************
 FUNCTION IO_hybridIA(Nast,ODFfileName)

 use prec, only: pReal, pInt
 use math, only: inRad

 implicit none
 
 character(len=*) ODFfileName
 character(len=80) line
 character(len=*), parameter :: fileFormat = '(A80)'
 integer(pInt) i,j,bin,Nast,NnonZero,Nset,Nreps,reps,phi1,Phi,phi2
 integer(pInt), dimension(7) :: pos
 integer(pInt), dimension(3) :: steps
 integer(pInt), dimension(:), allocatable :: binSet
 real(pReal) center,sum_dV_V,prob,dg_0,C,lowerC,upperC,rnd
 real(pReal), dimension(3) :: limits,deltas
 real(pReal), dimension(:,:,:), allocatable :: dV_V
 real(pReal), dimension(3,Nast) :: IO_hybridIA

 if (.not. IO_open_file(999,ODFfileName)) goto 100
 
!--- parse header of ODF file ---
!--- limits in phi1, Phi, phi2 ---
 read(999,fmt=fileFormat,end=100) line
 pos = IO_stringPos(line,3)
 if (pos(1).ne.3) goto 100
 do i=1,3
   limits(i) = IO_intValue(line,pos,i)*inRad
 end do

!--- deltas in phi1, Phi, phi2 ---
 read(999,fmt=fileFormat,end=100) line
 pos = IO_stringPos(line,3)
 if (pos(1).ne.3) goto 100
 do i=1,3
   deltas(i) = IO_intValue(line,pos,i)*inRad
 end do
 steps = nint(limits/deltas,pInt)
 allocate(dV_V(steps(3),steps(2),steps(1)))

!--- box boundary/center at origin? ---
 read(999,fmt=fileFormat,end=100) line
 if (index(IO_lc(line),'bound')>0) then
   center = 0.5_pReal
 else
   center = 0.0_pReal
 end if
 
!--- skip blank line ---
 read(999,fmt=fileFormat,end=100) line

 sum_dV_V = 0.0_pReal
 dV_V = 0.0_pReal
 dg_0 = deltas(1)*deltas(3)*2.0_pReal*sin(deltas(2)/2.0_pReal)
 NnonZero = 0_pInt
 
 do phi1=1,steps(1)
   do Phi=1,steps(2)
     do phi2=1,steps(3)
       read(999,fmt=*,end=100) prob
       if (prob > 0.0_pReal) then
         NnonZero = NnonZero+1
         sum_dV_V = sum_dV_V+prob
       else
         prob = 0.0_pReal
       end if
       dV_V(phi2,Phi,phi1) = prob*dg_0*sin((Phi-1.0_pReal+center)*deltas(2))
     end do
   end do
 end do  
 
 dV_V = dV_V/sum_dV_V  ! normalize to 1
 
!--- now fix bounds ---
 Nset = max(Nast,NnonZero)
 lowerC = 0.0_pReal
 upperC = real(Nset, pReal)
 
 do while (hybridIA_reps(dV_V,steps,upperC) < Nset)
   lowerC = upperC
   upperC = upperC*2.0_pReal
 end do
!--- binary search for best C ---
 do
   C = (upperC+lowerC)/2.0_pReal
   Nreps = hybridIA_reps(dV_V,steps,C)
   if (abs(upperC-lowerC) < upperC*1.0e-14_pReal) then
     C = upperC
     Nreps = hybridIA_reps(dV_V,steps,C)
     exit
   elseif (Nreps < Nset) then
     lowerC = C
   elseif (Nreps > Nset) then
     upperC = C
   else
     exit
   end if
 end do
 
 allocate(binSet(Nreps))
 bin = 0 ! bin counter
 i = 1 ! set counter
 do phi1=1,steps(1)
   do Phi=1,steps(2)
     do phi2=1,steps(3)
       reps = nint(C*dV_V(phi2,Phi,phi1), pInt)
       binSet(i:i+reps-1) = bin
       bin = bin+1 ! advance bin
       i = i+reps ! advance set
     end do
   end do
 end do

 do i=1,Nast
   if (i < Nast) then
     call random_number(rnd)
     j = nint(rnd*(Nast-i)+i+0.5_pReal,pInt)
   else
     j = i
   end if
   bin = binSet(j)
   IO_hybridIA(1,i) = deltas(1)*(mod(bin/(steps(3)*steps(2)),steps(1))+center)  ! phi1
   IO_hybridIA(2,i) = deltas(2)*(mod(bin/ steps(3)          ,steps(2))+center)  ! Phi
   IO_hybridIA(3,i) = deltas(3)*(mod(bin                    ,steps(3))+center)  ! phi2
   binSet(j) = binSet(i)
 end do
 close(999)
 return

! on error
100 IO_hybridIA = -1
 close(999)
 return
 
 END FUNCTION 


!********************************************************************
! identifies lines without content
!********************************************************************
 PURE FUNCTION IO_isBlank (line)

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
 
 END FUNCTION

!********************************************************************
! get tagged content of line
!********************************************************************
 PURE FUNCTION IO_getTag (line,openChar,closechar)

 use prec, only: pInt
 implicit none

 character(len=*), intent(in) :: line,openChar,closeChar
 character(len=*), parameter :: sep=achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces
 character(len=len_trim(line)) IO_getTag
 integer(pInt)  left,right

 IO_getTag = ''
 left = scan(line,openChar)
 right = scan(line,closeChar)
 
 if (left == verify(line,sep) .and. right > left) & ! openChar is first and closeChar occurs
   IO_getTag = line(left+1:right-1)

 return
 
 END FUNCTION


!*********************************************************************
 FUNCTION IO_countSections(file,part)
!*********************************************************************
 use prec, only: pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: file
 character(len=*), intent(in) :: part
 integer(pInt) IO_countSections
 character(len=1024) line

 IO_countSections = 0
 line = ''
 rewind(file)

 do while (IO_getTag(line,'<','>') /= part)      ! search for part
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     IO_countSections = IO_countSections + 1
 enddo

100 return

 END FUNCTION
 

!*********************************************************************
! return array of myTag counts within <part> for at most N[sections]
!*********************************************************************
 FUNCTION IO_countTagInPart(file,part,myTag,Nsections)

 use prec, only: pInt
 implicit none

!* Definition of variables
 integer(pInt), intent(in) :: file, Nsections
 character(len=*), intent(in) :: part, myTag
 integer(pInt), dimension(Nsections) :: IO_countTagInPart, counter
 integer(pInt), parameter :: maxNchunks = 1
 integer(pInt), dimension(1+2*maxNchunks) :: positions
 integer(pInt) section
 character(len=1024) line,tag

 counter = 0_pInt
 section = 0_pInt
 line = ''
 rewind(file)

 do while (IO_getTag(line,'<','>') /= part)               ! search for part
   read(file,'(a1024)',END=100) line
 enddo

 do
   read(file,'(a1024)',END=100) line
   if (IO_isBlank(line)) cycle                            ! skip empty lines
   if (IO_getTag(line,'<','>') /= '') exit                ! stop at next part
   if (IO_getTag(line,'[',']') /= '') &                   ! found [section] identifier
     section = section + 1
   if (section > 0) then
     positions = IO_stringPos(line,maxNchunks)
     tag = IO_lc(IO_stringValue(line,positions,1))        ! extract key
     if (tag == myTag) &                                  ! match
       counter(section) = counter(section) + 1
   endif   
 enddo

100 IO_countTagInPart = counter
 return

END FUNCTION


!********************************************************************
! locate at most N space-separated parts in line
! return array containing number of parts found and
! their left/right positions to be used by IO_xxxVal
!********************************************************************
 PURE FUNCTION IO_stringPos (line,N)

 use prec, only: pReal,pInt
 implicit none

 character(len=*), intent(in) :: line
 character(len=*), parameter :: sep=achar(32)//achar(9)//achar(10)//achar(13) ! whitespaces
 integer(pInt), intent(in) :: N
 integer(pInt)  part
 integer(pInt) IO_stringPos(1+N*2)

 IO_stringPos = -1
 IO_stringPos(1) = 0
 part = 1
 do while ((N<1 .or. part<=N) .and. verify(line(IO_stringPos(part*2-1)+1:),sep)>0)
   IO_stringPos(part*2) = IO_stringPos(part*2-1)+verify(line(IO_stringPos(part*2-1)+1:),sep)
   IO_stringPos(part*2+1) = IO_stringPos(part*2)+scan(line(IO_stringPos(part*2):),sep)-2
   part = part+1
 end do
 IO_stringPos(1) = part-1
 return

 END FUNCTION


!********************************************************************
! read string value at pos from line
!********************************************************************
 PURE FUNCTION IO_stringValue (line,positions,pos)
 
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

 END FUNCTION


!********************************************************************
! read string value at pos from fixed format line
!********************************************************************
 PURE FUNCTION IO_fixedStringValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 character(len=ends(pos+1)-ends(pos)) IO_fixedStringValue

 IO_fixedStringValue = line(ends(pos)+1:ends(pos+1))
 return

 END FUNCTION


!********************************************************************
! read float value at pos from line
!********************************************************************
 PURE FUNCTION IO_floatValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 real(pReal) IO_floatValue

 if (positions(1) >= pos) then
   read(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT=*) IO_floatValue
   return
 endif
100 IO_floatValue = huge(1.0_pReal)
 return

 END FUNCTION


!********************************************************************
! read float value at pos from fixed format line
!********************************************************************
 PURE FUNCTION IO_fixedFloatValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 real(pReal) IO_fixedFloatValue

 read(UNIT=line(ends(pos-1)+1:ends(pos)),ERR=100,FMT=*) IO_fixedFloatValue
 return
100 IO_fixedFloatValue = huge(1.0_pReal)
 return

 END FUNCTION


!********************************************************************
! read float x.y+z value at pos from format line line
!********************************************************************
 PURE FUNCTION IO_fixedNoEFloatValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 integer(pInt) pos_exp,expon
 real(pReal) IO_fixedNoEFloatValue,base
 
 pos_exp = scan(line(ends(pos)+1:ends(pos+1)),'+-',back=.true.)
 if (pos_exp > 1) then
   read(UNIT=line(ends(pos)+1:ends(pos)+pos_exp-1),ERR=100,FMT=*) base
   read(UNIT=line(ends(pos)+pos_exp:ends(pos+1)),ERR=100,FMT=*) expon
 else
   read(UNIT=line(ends(pos)+1:ends(pos+1)),ERR=100,FMT=*) base
   expon = 0_pInt
 endif
 IO_fixedNoEFloatValue = base*10.0_pReal**expon
 return
100 IO_fixedNoEFloatValue = huge(1.0_pReal)
 return

 END FUNCTION


!********************************************************************
! read int value at pos from line
!********************************************************************
 PURE FUNCTION IO_intValue (line,positions,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: positions(*),pos
 integer(pInt) IO_intValue

 if (positions(1) >= pos) then
   read(UNIT=line(positions(pos*2):positions(pos*2+1)),ERR=100,FMT=*) IO_intValue
   return
 endif
100 IO_intValue = huge(1_pInt)
 return

 END FUNCTION


!********************************************************************
! read int value at pos from fixed format line
!********************************************************************
 PURE FUNCTION IO_fixedIntValue (line,ends,pos)
 
 use prec, only: pReal,pInt
 implicit none
 
 character(len=*), intent(in) :: line
 integer(pInt), intent(in) :: ends(*),pos
 integer(pInt) IO_fixedIntValue

 read(UNIT=line(ends(pos)+1:ends(pos+1)),ERR=100,FMT=*) IO_fixedIntValue
 return
100 IO_fixedIntValue = huge(1_pInt)
 return

 END FUNCTION


!********************************************************************
! change character in line to lower case
!********************************************************************
 PURE FUNCTION IO_lc (line)

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

 END FUNCTION


!********************************************************************
! in place change of character in line to lower case
!********************************************************************
 SUBROUTINE IO_lcInplace (line)

 use prec, only: pInt
 implicit none

 character (len=*) line
 character (len=len(line)) IO_lc
 integer(pInt) i

 IO_lc = line
 do i=1,len(line)
    if(64<iachar(line(i:i)) .and. iachar(line(i:i))<91) IO_lc(i:i)=achar(iachar(line(i:i))+32)
 enddo
 line = IO_lc
 return 

 END SUBROUTINE


!********************************************************************
! count items in consecutive lines of ints concatenated by "c"
! as last char or range of values a "to" b
!********************************************************************
 FUNCTION IO_countContinousIntValues (unit)

 use prec, only: pReal,pInt
 implicit none

 integer(pInt)  IO_countContinousIntValues,unit
 integer(pInt), dimension(67) :: pos  ! allow for 32 values excl "c"
 character(len=300) line

 IO_countContinousIntValues = 0
 do
   read(unit,'(A300)',end=100) line
   pos = IO_stringPos(line,33)
   if (IO_lc(IO_stringValue(line,pos,2)) == 'to' ) then  ! found range indicator
     IO_countContinousIntValues = IO_countContinousIntValues+1+IO_intValue(line,pos,3)-IO_intValue(line,pos,1)
	 exit
   else
     IO_countContinousIntValues = IO_countContinousIntValues+pos(1)-1
     if ( IO_lc(IO_stringValue(line,pos,pos(1))) /= 'c' ) then  ! line finished, read last value
       IO_countContinousIntValues = IO_countContinousIntValues+1
       exit
     endif
   endif
 enddo
100 return

 END FUNCTION

!*********************************************************************
! read consecutive lines of ints concatenated by "c" as last char
! or range of values a "to" b
!*********************************************************************
 FUNCTION IO_continousIntValues (unit,maxN,lookupName,lookupMap,lookupMaxN)

 use prec, only: pReal,pInt
 implicit none

 integer(pInt)  unit,maxN,i
 integer(pInt), dimension(1+maxN) :: IO_continousIntValues
 integer(pInt), dimension(67) :: pos  ! allow for 32 values excl "c"
 character(len=64), dimension(:) :: lookupName
 integer(pInt) :: lookupMaxN
 integer(pInt), dimension(:,:) :: lookupMap
 character(len=300) line

 IO_continousIntValues = 0_pInt
 do
   read(unit,'(A300)',end=100) line
   pos = IO_stringPos(line,33)
   if (verify(IO_stringValue(line,pos,1),"0123456789") > 0) then     ! a non-int, i.e. set name
     do i = 1,lookupMaxN                                       ! loop over known set names
       if (IO_stringValue(line,pos,1) == lookupName(i)) then   ! found matching name
         IO_continousIntValues = lookupMap(:,i)                ! return resp. entity list
         exit
       endif
     enddo
     exit
   else if (IO_lc(IO_stringValue(line,pos,2)) == 'to' ) then   ! found range indicator
     do i = IO_intValue(line,pos,1),IO_intValue(line,pos,3)
       IO_continousIntValues(1) = IO_continousIntValues(1)+1
	   IO_continousIntValues(1+IO_continousIntValues(1)) = i
     enddo
     exit
   else
     do i = 1,pos(1)-1  ! interpret up to second to last value
       IO_continousIntValues(1) = IO_continousIntValues(1)+1
       IO_continousIntValues(1+IO_continousIntValues(1)) = IO_intValue(line,pos,i)
     enddo
     if ( IO_lc(IO_stringValue(line,pos,pos(1))) /= 'c' ) then  ! line finished, read last value
       IO_continousIntValues(1) = IO_continousIntValues(1)+1
       IO_continousIntValues(1+IO_continousIntValues(1)) = IO_intValue(line,pos,pos(1))
       exit
     endif
   endif
 enddo
100 return

 END FUNCTION


!********************************************************************
! write error statements to standard out
! and terminate the Marc run with exit #9xxx
! in ABAQUS either time step is reduced or execution terminated
!********************************************************************
 SUBROUTINE IO_error(ID,e,i,g,ext_msg)

 use prec, only: pInt

 use debug
 implicit none

 integer(pInt), intent(in) :: ID
 integer(pInt), optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 character(len=80) msg

 select case (ID)
 case (0)
   msg = 'Unable to open input file'
 case (100)
   msg = 'Error reading from configuration file'
 case (105)
   msg = 'Error reading from ODF file'
 case (110)
   msg = 'No homogenization specified via State Variable 2'
 case (120)
   msg = 'No microstructure specified via State Variable 3'
 case (130)
   msg = 'Homogenization index out of bounds'
 case (140)
   msg = 'Microstructure index out of bounds'
 case (150)
   msg = 'Phase index out of bounds'
 case (160)
   msg = 'Texture index out of bounds'
 case (170)
   msg = 'Sum of phase fractions differs from 1'
 case (200)
   msg = 'Unknown constitution specified'
 case (201)
   msg = 'Unknown lattice type specified'
 case (202)
   msg = 'Number of slip systems too small'
 case (203)
   msg = 'Negative initial slip resistance'
 case (204)
   msg = 'Non-positive reference shear rate'
 case (205)
   msg = 'Non-positive stress exponent'
 case (206)
   msg = 'Non-positive initial hardening slope'
 case (207)
   msg = 'Non-positive saturation stress'
 case (208)
   msg = 'Non-positive w0'
 case (209)
   msg = 'Negative latent hardening ratio'
 case (220)
   msg = 'Negative initial dislocation density'
 case (221)
   msg = 'Negative Bugers vector'
 case (222)
   msg = 'Negative activation energy for edge dislocation glide'
 case (223)
   msg = 'Negative self diffusion energy'
 case (224)
   msg = 'Negative diffusion constant'
 case (240)
   msg = 'Non-positive Taylor factor'
 case (300)
   msg = 'This material can only be used with elements with three direct stress components'
 case (500)
   msg = 'Unknown lattice type specified'
 case (600)
   msg = 'Convergence not reached'
 case (610)
   msg = 'Stress loop not converged'
 case (700)
   msg = 'Singular matrix in stress iteration'
 case (800)
   msg = 'GIA requires 8 grains per IP (bonehead, you!)'
 case default
   msg = 'Unknown error number...'
 end select
 
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '+------------------------------+'
 write(6,*) '+            ERROR             +'
 write(6,*) '+                              +'
 write(6,*) msg
 if (present(ext_msg)) write(6,*) ext_msg
 if (present(e)) then
   if (present(i) .and. present(g)) then
     write(6,'(a10,x,i6,x,a2,x,i2,x,a5,x,i4)') 'at element',e,'IP',i,'grain',g
   else
     write(6,'(a2,x,i6)') 'at',e
   endif
 endif
 write(6,*) '+------------------------------+'

 call debug_info()
 call flush(6)
 call quit(9000+ID)
!$OMP END CRITICAL (write2out)

! ABAQUS returns in some cases
 return

 END SUBROUTINE


!********************************************************************
! write warning statements to standard out
!********************************************************************
 SUBROUTINE IO_warning(ID,e,i,g,ext_msg)

 use prec, only: pInt
 use debug
 implicit none

 integer(pInt), intent(in) :: ID
 integer(pInt), optional, intent(in) :: e,i,g
 character(len=*), optional, intent(in) :: ext_msg
 character(len=80) msg

 select case (ID)
 case (650)
   msg = 'Polar decomposition failed'
 case default
   msg = 'Unknown warning number...'
 end select
 
 !$OMP CRITICAL (write2out)
 write(6,*)
 write(6,*) '+------------------------------+'
 write(6,*) '+           warning            +'
 write(6,*) '+                              +'
 write(6,*) msg
 if (present(ext_msg))  write(6,*) ext_msg
 if (present(e)) then
   if (present(i) .and. present(g)) then
     write(6,'(a10,x,i6,x,a2,x,i2,x,a5,x,i4)') 'at element',e,'IP',i,'grain',g
   else
     write(6,'(a2,x,i6)') 'at',e
   endif
 endif
 write(6,*) '+------------------------------+'

 END SUBROUTINE
 
 END MODULE IO
