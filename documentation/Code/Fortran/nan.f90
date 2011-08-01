 program nantest
! This file: http://ftp.aset.psu.edu/pub/ger/fortran/hdk/nan.f90
!
!====NAN.F90 illustrates what works and what doesn't when
!    detecting a NaN. All but the Fortran 2003 IEEE_IS_NAN are not
!    Standard Fortran.
!
! Platforms: Windows 9x/Me/NT/2000, AIX 4.3, Linux
! Compiler Notes:
!  Compaq Visual Fortran 6.6C with default
!     fpe settings (/fpe:3 /nocheck) and /OPTIMIZE:0
!     (ISNAN is an Elemental Intrinsic Function)
!     Options /fpe:0 /traceback will cause this
!     program to stop with error message,
!     James Giles points out that in order to actually print
!     minus zero from CVF, you have to compile with the
!     /assume:minus0  option.
!
!   AIX XLF90 without optimization.
!   (ISNAN is part of a BOS Runtime C library;
!   thus ISNAN must be declared LOGICAL.)
!
!   Linux Intel IFORT, use: -O0 (no optimization) and
!                           -assume minus0
!
!   g95 supports ieee_is_nan elemental function.
!
!
! Author: hdkLESS at SPAM psu dot edu
! Date: March, 2002, January 2005, August 2006.
!
! References:
! http://www.psc.edu/general/software/packages/ieee/ieee.html
! http://homepages.borland.com/efg2lab/Mathematics/NaN.htm
! http://en.wikipedia.org/wiki/NaN
!
       logical :: ISNAN
       integer :: i
       integer, Dimension (6) :: Iy
       real, Dimension (6) :: y
       integer :: IPInf, IMinf, IMZero
       real ::  PInf, MInf, MZero, DivNan
       Character (Len=10), Dimension(6) :: NType
       data NType/'+Infinity','-Infinity','-0', 'NaN','NaN','0/0'/
       data IPInf/B'01111111100000000000000000000000'/    ! +Infinity
       data IMInf/B'11111111100000000000000000000000'/    ! -Infinity
       data IMZero/B'10000000000000000000000000000000'/   ! -0

       data Iy(1)/B'01111111100000000000000000000000'/       ! +Infinity
       data Iy(2)/B'11111111100000000000000000000000'/       ! -Infinity
       data Iy(3)/B'10000000000000000000000000000000'/       ! -0
       data Iy(4)/B'01111111100000100000000000000000'/       ! NaN
       data Iy(5)/B'11111111100100010001001010101010'/       ! NaN
       data Iy(6)/B'11111111110000000000000000000000'/       ! 0/0

       PInf = transfer(IPinf,Pinf)
       Minf = transfer(IMinf,Minf)
       MZero = transfer(IMZero,MZero)
       Y = transfer(IY,Y)

! Some debug options for some compilers may flag the following
! zero/zero as an exception. If so, comment out the next two lines.
        DivNan=0
        y(6)=DivNan/DivNan

       Do i=1,6
        print *, 'Test#',i,' ->',NType(i)
        if (y(i).eq.PInf) print *, 'Y = Plus Infinity'
        if (y(i).eq.MInf) print *, 'Y = Minus Infinity'
        if (y(i).eq.Mzero) print *, 'Y = Minus Zero'
        print *, 'y(Test#)=',y(i)
        print *, 'Testing each of three NaN detection methods:'
! EQV -> true iff both A and B are True or iff both A and B are False.
        if( (y(i) > 0.0) .EQV. (y(i) <= 0.0)) then
           print *, '1) (y(Test#) > 0.0) .EQV. (y(Test#) <= 0.0)'
        end if
        if (y(i)/=y(i)) then
           print *, '2) (y(Test#)/=(y(Test#))'
        end if
! If ISNAN is available for a specific compiler, uncomment the
! following 3 lines.
!         if (ISNAN(y(i))) then
!           print *, '3) ISNAN(y(Test#))'
!         end if
! if Fortran 2003 IEEE floating-point exception handling is available
! uncomment the following 4 lines
!        use ieee_arithmetic
!        if (ieee_is_nan(y(i))) then
!          print *, '4) ieee_is_nan(y(Test#))'
!        end if
         print *, ' '
       End Do

 end program nantest
