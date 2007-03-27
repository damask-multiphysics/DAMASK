! MISSING agree on Rad as standard in/out
! MISSING inverse 3x3 and 6x6

!##############################################################
 MODULE math   
!##############################################################

 use prec, only: pReal,pInt
 implicit none

 real(pReal), parameter :: pi = 3.14159265358979323846264338327950288419716939937510_pReal
 real(pReal), parameter :: inDeg = 180.0_pReal/pi
 real(pReal), parameter :: inRad = pi/180.0_pReal
 real(pReal), dimension(3,3), parameter :: math_I3 = &
 reshape( (/ &
 1.0_pReal,0.0_pReal,0.0_pReal, &
 0.0_pReal,1.0_pReal,0.0_pReal, &
 0.0_pReal,0.0_pReal,1.0_pReal /),(/3,3/))

 CONTAINS

 SUBROUTINE math_invert3x3(A, InvA, DetA, err)

!   Bestimmung der Determinanten und Inversen einer 3x3-Matrix

!   A      = Matrix A
!   InvA   = Inverse von A
!   DetA   = Determinante von A
!   Fehler = Fehlerparameter
!            0 : Die Berechnung wurde durchgefuehrt.
!            1 : Die Determinante ist kleiner gleich Null.

 use prec, only: pReal,pInt
 implicit none

 integer (pInt) err
 real(pReal)    InvA(3,3), DetA, A(3,3)

 INTENT (IN)     A
 INTENT (OUT)    InvA, DetA, err

 err = 0

 DetA =   A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) )&
        - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) )&
        + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) )

 IF (DetA <= 0.0000001) THEN
 err = 1
 RETURN
 END IF

 InvA(1,1) = (  A(2,2) * A(3,3) - A(2,3) * A(3,2) ) / DetA
 InvA(2,1) = ( -A(2,1) * A(3,3) + A(2,3) * A(3,1) ) / DetA
 InvA(3,1) = (  A(2,1) * A(3,2) - A(2,2) * A(3,1) ) / DetA

 InvA(1,2) = ( -A(1,2) * A(3,3) + A(1,3) * A(3,2) ) / DetA
 InvA(2,2) = (  A(1,1) * A(3,3) - A(1,3) * A(3,1) ) / DetA
 InvA(3,2) = ( -A(1,1) * A(3,2) + A(1,2) * A(3,1) ) / DetA

 InvA(1,3) = (  A(1,2) * A(2,3) - A(1,3) * A(2,2) ) / DetA
 InvA(2,3) = ( -A(1,1) * A(2,3) + A(1,3) * A(2,1) ) / DetA
 InvA(3,3) = (  A(1,1) * A(2,2) - A(1,2) * A(2,1) ) / DetA

 RETURN

 END SUBROUTINE math_invert3x3
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE math_invert6x6(A, InvA, AnzNegEW, err)

!   Invertieren einer 6x6-Matrix

!   A        = Matrix A
!   InvA     = Inverse von A
!   AnzNegEW = Anzahl der negativen Eigenwerte von A
!   Fehler   = Fehlerparameter
!              = 0: Inversion wurde durchgefuehrt.
!              = 1: Die Inversion in SymGauss wurde wegen eines verschwindenen
!                   Pivotelement abgebrochen.

 use prec, only: pReal,pInt
 implicit none

 integer (pInt) AnzNegEW, err
 real(pReal)    InvA(6,6), A(6,6)


 INTENT (IN)     A
 INTENT (OUT)    InvA, AnzNegEW, err


 integer (pInt) i
 real(pReal)    LogAbsDetA
 real(pReal)    B(6,6)

 InvA = 0.0_pReal

 forall(i = 1:6) InvA(i,i) = 1.0_pReal
 B = A
 CALL Gauss (B, InvA, LogAbsDetA, AnzNegEW, err)

 RETURN

 END SUBROUTINE math_invert6x6

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

      SUBROUTINE Gauss (A, B,LogAbsDetA, NegHDK,Fehler)

!   Loesung eines linearen Gleichungsssystem A * X = B mit Hilfe des
!   GAUSS-Algorithmusses

!   Zur numerischen Stabilisierung wird eine Zeilen- und Spaltenpivotsuche
!   durchgefuehrt.

!   Eingabeparameter:
!
!   A(6,6) = Koeffizientenmatrix A
!   B(6,6) = rechte Seiten B
!
!   Ausgabeparameter:
!
!   B(6,6)        = Matrix der Unbekanntenvektoren X
!   LogAbsDetA    = 10-Logarithmus des Betrages der Determinanten von A
!   NegHDK        = Anzahl der negativen Hauptdiagonalkoeffizienten nach der
!                   Vorwaertszerlegung
!   Fehler        = Fehlerparameter
!                   = 0: Das Gleichungssystem wurde geloest.
!                   = 1: Matrix A ist singulaer.

!   A und B werden veraendert!

 use prec, only: pReal,pInt
 implicit none

 integer (pInt) NegHDK, Fehler
 real(pReal)    LogAbsDetA
 real(pReal)    A(6,6), B(6,6)

 INTENT (OUT)    LogAbsDetA, NegHDK, Fehler
 INTENT (INOUT)  A, B

 LOGICAL        SortX
 integer (pInt) PivotZeile, PivotSpalte, StoreI, I, IP1, J, K, L
 integer (pInt) XNr(6)
 real(pReal)    AbsA, PivotWert, EpsAbs, Quote
 real(pReal)    StoreA(6), StoreB(6)

 Fehler = 1
 NegHDK = 1
 SortX  = .FALSE.

!   Unbekanntennumerierung

 DO  I = 1, 6
    XNr(I) = I
 END DO

!   Genauigkeitsschranke und Bestimmung des groessten Pivotelementes

 PivotWert   = ABS(A(1,1))
 PivotZeile  = 1
 PivotSpalte = 1

 DO  I = 1, 6
    DO  J = 1, 6
        AbsA = ABS(A(I,J))
        IF (AbsA .GT. PivotWert) THEN
            PivotWert   = AbsA
            PivotZeile  = I
            PivotSpalte = J
        END IF
    END DO
 END DO

 IF (PivotWert .LT. 0.0000001) RETURN   ! Pivotelement = 0?

 EpsAbs = PivotWert * 0.1_pReal ** PRECISION(1.0_pReal)

!   V O R W A E R T S T R I A N G U L A T I O N

 DO  I = 1, 6 - 1

!     Zeilentausch?

    IF (PivotZeile .NE. I) THEN

        StoreA(I:6)       = A(I,I:6)
        A(I,I:6)          = A(PivotZeile,I:6)
        A(PivotZeile,I:6) = StoreA(I:6)

        StoreB(1:6)        = B(I,1:6)
        B(I,1:6)           = B(PivotZeile,1:6)
        B(PivotZeile,1:6)  = StoreB(1:6)

        SortX                = .TRUE.

    END IF

!     Spaltentausch?

    IF (PivotSpalte .NE. I) THEN

        StoreA(1:6)        = A(1:6,I)
        A(1:6,I)           = A(1:6,PivotSpalte)
        A(1:6,PivotSpalte) = StoreA(1:6)

        StoreI                = XNr(I)
        XNr(I)                = XNr(PivotSpalte)
        XNr(PivotSpalte)      = StoreI

        SortX                 = .TRUE.

    END IF

!     Triangulation

    DO  J = I + 1, 6

        Quote = A(J,I) / A(I,I)

        DO  K = I + 1, 6
            A(J,K) = A(J,K) - Quote * A(I,K)
        END DO

        DO  K = 1, 6
            B(J,K) = B(J,K) - Quote * B(I,K)
        END DO

    END DO

!     Bestimmung des groessten Pivotelementes

    IP1         = I + 1
    PivotWert   = ABS(A(IP1,IP1))
    PivotZeile  = IP1
    PivotSpalte = IP1

    DO  J = IP1, 6
        DO  K = IP1, 6
            AbsA = ABS(A(J,K))
            IF (AbsA .GT. PivotWert) THEN
                PivotWert   = AbsA
                PivotZeile  = J
                PivotSpalte = K
            END IF
        END DO
    END DO

    IF (PivotWert .LT. EpsAbs) RETURN   ! Pivotelement = 0?

 END DO

!   R U E C K W A E R T S A U F L O E S U N G

 DO  I = 6, 1, -1

    DO  L = 1, 6

        DO  J = I + 1, 6
            B(I,L) = B(I,L) - A(I,J) * B(J,L)
        END DO

        B(I,L) = B(I,L) / A(I,I)

    END DO

 END DO

!   Sortieren der Unbekanntenvektoren?

 IF (SortX) THEN

    DO  L = 1, 6

        StoreA(1:6) = B(1:6,L)

        DO  I = 1, 6
            J      = XNr(I)
            B(J,L) = StoreA(I)
        END DO

    END DO

 END IF

!   Determinante

 LogAbsDetA = 0.0_pReal
 NegHDK     = 0

 DO  I = 1, 6
    IF (A(I,I) .LT. 0.0_pReal) NegHDK = NegHDK + 1
    AbsA       = ABS(A(I,I))
    LogAbsDetA = LogAbsDetA + LOG10(AbsA)
 END DO

 Fehler = 0

 RETURN

 END SUBROUTINE Gauss





!    ***  Initialize random number generator     ***
!    *** for later use in mpie_fiber and mpie_disturbOri ***

 SUBROUTINE math_init ()

 use prec, only: pReal,pInt
 implicit none

 integer (pInt) seed
 
 call random_seed()
 call get_seed(seed)
 call halton_seed_set(seed)
 call halton_ndim_set(3)

 END SUBROUTINE
 
 
!********************************************************************
! calculate v6 stress from M6x6 stiffness and v6 strain
!********************************************************************
 FUNCTION math_Hooke(C,e)

 use prec, only: pReal,pInt
 implicit none

 real(pReal) C(6,6), e(6), factoredE(6)
 real(pReal), dimension(6), parameter :: preFactor = (/1.0_pReal,1.0_pReal,1.0_pReal,2.0_pReal,2.0_pReal,2.0_pReal/)
 real(pReal) math_Hooke(6)

 factoredE = dot_product(preFactor,e)
 math_Hooke = matmul(C,factoredE)
 return

 END FUNCTION

 
!********************************************************************
! calculate the determinant of a (3x3)
!********************************************************************
 real(pReal) FUNCTION math_det3x3(m)

 use prec, only: pReal,pInt
 implicit none

 real(pReal) m(3,3)

 math_det3x3 = m(1,1)*(m(2,2)*m(3,3)-m(2,3)*m(3,2)) &
              -m(1,2)*(m(2,1)*m(3,3)-m(2,3)*m(3,1)) &
              +m(1,3)*(m(2,1)*m(3,2)-m(2,2)*m(3,1))
 return

 END FUNCTION


!********************************************************************
! convert a symmetric 3,3 matrix into an array of 6
!********************************************************************
 FUNCTION math_33to6(m33)

 use prec, only: pReal,pInt
 implicit none

 real(pReal), dimension(3,3) :: m33
 real(pReal), dimension(6) :: math_33to6

 math_33to6(1)=m33(1,1)
 math_33to6(2)=m33(2,2)
 math_33to6(3)=m33(3,3)
 math_33to6(4)=m33(1,2)
 math_33to6(5)=m33(2,3)
 math_33to6(6)=m33(1,3)
 return

 END FUNCTION


!********************************************************************
! This routine coverts an array of 6 into a symmetric 3,3 matrix
!********************************************************************
 FUNCTION math_6to33(v6)

 use prec, only: pReal,pInt
 implicit none

 real(pReal) math_6to33(3,3), v6(6)

 math_6to33(1,1)=v6(1)
 math_6to33(2,2)=v6(2)
 math_6to33(3,3)=v6(3)
 math_6to33(1,2)=v6(4)
 math_6to33(2,1)=v6(4)
 math_6to33(2,3)=v6(5)
 math_6to33(3,2)=v6(5)
 math_6to33(1,3)=v6(6)
 math_6to33(3,1)=v6(6)
 return

 END FUNCTION

 
!********************************************************************************
!**      This routine transforms the stiffness matrix    **
!********************************************************************************
 FUNCTION math_66to3333(C66)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) C66(6,6), math_66to3333(3,3,3,3)

 math_66to3333(1,1,1,1)=C66(1,1)
 math_66to3333(1,1,2,2)=C66(1,2)
 math_66to3333(1,1,3,3)=C66(1,3) 
 math_66to3333(1,1,2,3)=C66(1,4)
 math_66to3333(1,1,3,2)=C66(1,4)
 math_66to3333(1,1,1,3)=C66(1,5) 
 math_66to3333(1,1,3,1)=C66(1,5)
 math_66to3333(1,1,1,2)=C66(1,6)
 math_66to3333(1,1,2,1)=C66(1,6) 
 math_66to3333(2,2,1,1)=C66(2,1)
 math_66to3333(2,2,2,2)=C66(2,2)
 math_66to3333(2,2,3,3)=C66(2,3) 
 math_66to3333(2,2,2,3)=C66(2,4)
 math_66to3333(2,2,3,2)=C66(2,4)
 math_66to3333(2,2,1,3)=C66(2,5) 
 math_66to3333(2,2,3,1)=C66(2,5)
 math_66to3333(2,2,1,2)=C66(2,6)
 math_66to3333(2,2,2,1)=C66(2,6) 
 math_66to3333(3,3,1,1)=C66(3,1)
 math_66to3333(3,3,2,2)=C66(3,2)
 math_66to3333(3,3,3,3)=C66(3,3) 
 math_66to3333(3,3,2,3)=C66(3,4)
 math_66to3333(3,3,3,2)=C66(3,4)
 math_66to3333(3,3,1,3)=C66(3,5) 
 math_66to3333(3,3,3,1)=C66(3,5)
 math_66to3333(3,3,1,2)=C66(3,6)
 math_66to3333(3,3,2,1)=C66(3,6) 
 math_66to3333(2,3,1,1)=C66(4,1)
 math_66to3333(3,2,1,1)=C66(4,1)
 math_66to3333(2,3,2,2)=C66(4,2) 
 math_66to3333(3,2,2,2)=C66(4,2)
 math_66to3333(2,3,3,3)=C66(4,3)
 math_66to3333(3,2,3,3)=C66(4,3)
 math_66to3333(2,3,2,3)=C66(4,4)
 math_66to3333(2,3,3,2)=C66(4,4)
 math_66to3333(3,2,2,3)=C66(4,4) 
 math_66to3333(3,2,3,2)=C66(4,4)
 math_66to3333(2,3,3,1)=C66(4,5)
 math_66to3333(2,3,1,3)=C66(4,5)
 math_66to3333(2,3,3,1)=C66(4,5)
 math_66to3333(3,2,1,3)=C66(4,5)
 math_66to3333(2,3,1,2)=C66(4,6) 
 math_66to3333(2,3,2,1)=C66(4,6)
 math_66to3333(3,2,1,2)=C66(4,6)
 math_66to3333(3,2,2,1)=C66(4,6) 
 math_66to3333(3,1,1,1)=C66(5,1)
 math_66to3333(1,3,1,1)=C66(5,1)
 math_66to3333(3,1,2,2)=C66(5,2) 
 math_66to3333(1,3,2,2)=C66(5,2)
 math_66to3333(3,1,3,3)=C66(5,3)
 math_66to3333(1,3,3,3)=C66(5,3)
 math_66to3333(3,1,2,3)=C66(5,4)
 math_66to3333(3,1,3,2)=C66(5,4)
 math_66to3333(1,3,2,3)=C66(5,4) 
 math_66to3333(1,3,3,2)=C66(5,4)
 math_66to3333(3,1,3,1)=C66(5,5)
 math_66to3333(3,1,1,3)=C66(5,5) 
 math_66to3333(1,3,3,1)=C66(5,5)
 math_66to3333(1,3,1,3)=C66(5,5)
 math_66to3333(3,1,1,2)=C66(5,6)
 math_66to3333(3,1,2,1)=C66(5,6)
 math_66to3333(1,3,1,2)=C66(5,6)
 math_66to3333(1,3,2,1)=C66(5,6)
 math_66to3333(1,2,1,1)=C66(6,1)
 math_66to3333(2,1,1,1)=C66(6,1)
 math_66to3333(1,2,2,2)=C66(6,2) 
 math_66to3333(2,1,2,2)=C66(6,2)
 math_66to3333(1,2,3,3)=C66(6,3)
 math_66to3333(2,1,3,3)=C66(6,3) 
 math_66to3333(1,2,2,3)=C66(6,4)
 math_66to3333(1,2,3,2)=C66(6,4)
 math_66to3333(2,1,2,3)=C66(6,4)
 math_66to3333(2,1,3,2)=C66(6,4)
 math_66to3333(1,2,3,1)=C66(6,5)
 math_66to3333(1,2,1,3)=C66(6,5) 
 math_66to3333(2,1,3,1)=C66(6,5)
 math_66to3333(2,1,1,3)=C66(6,5)
 math_66to3333(1,2,1,2)=C66(6,6)  
 math_66to3333(1,2,2,1)=C66(6,6)
 math_66to3333(2,1,1,2)=C66(6,6)
 math_66to3333(2,1,2,1)=C66(6,6)    
 return

 END FUNCTION 

     
 FUNCTION math_3333to66(C3333)
!********************************************************************************
!**      This routine transforms the stiffness matrix    **
!********************************************************************************
 use prec, only: pReal, pInt
 implicit none

 real(pReal) math_3333to66(6,6), C3333(3,3,3,3)

 math_3333to66(1,1)=C3333(1,1,1,1)
 math_3333to66(1,2)=C3333(1,1,2,2)
 math_3333to66(1,3)=C3333(1,1,3,3)
 math_3333to66(1,4)=C3333(1,1,2,3)
 math_3333to66(1,5)=C3333(1,1,3,1)
 math_3333to66(1,6)=C3333(1,1,1,2)
 math_3333to66(2,1)=C3333(2,2,1,1)
 math_3333to66(2,2)=C3333(2,2,2,2)
 math_3333to66(2,3)=C3333(2,2,3,3)
 math_3333to66(2,4)=C3333(2,2,2,3)
 math_3333to66(2,5)=C3333(2,2,3,1)
 math_3333to66(2,6)=C3333(2,2,1,2)
 math_3333to66(3,1)=C3333(3,3,1,1)
 math_3333to66(3,2)=C3333(3,3,2,2)
 math_3333to66(3,3)=C3333(3,3,3,3)
 math_3333to66(3,4)=C3333(3,3,2,3)
 math_3333to66(3,5)=C3333(3,3,3,1)
 math_3333to66(3,6)=C3333(3,3,1,2)
 math_3333to66(4,1)=C3333(2,3,1,1)
 math_3333to66(4,2)=C3333(2,3,2,2)
 math_3333to66(4,3)=C3333(2,3,3,3)
 math_3333to66(4,4)=C3333(2,3,2,3)
 math_3333to66(4,5)=C3333(2,3,3,1)
 math_3333to66(4,6)=C3333(2,3,1,2)
 math_3333to66(5,1)=C3333(3,1,1,1)
 math_3333to66(5,2)=C3333(3,1,2,2)
 math_3333to66(5,3)=C3333(3,1,3,3)
 math_3333to66(5,4)=C3333(3,1,2,3)
 math_3333to66(5,5)=C3333(3,1,3,1)
 math_3333to66(5,6)=C3333(3,1,1,2)  
 math_3333to66(6,1)=C3333(1,2,1,1)
 math_3333to66(6,2)=C3333(1,2,2,2)
 math_3333to66(6,3)=C3333(1,2,3,3) 
 math_3333to66(6,4)=C3333(1,2,2,3)
 math_3333to66(6,5)=C3333(1,2,3,1)
 math_3333to66(6,6)=C3333(1,2,1,2) 
 return

 END FUNCTION


!********************************************************************
! This routine calculates Euler angles from orientation matrix
!********************************************************************
 SUBROUTINE math_RtoEuler(orimat, phi1, PHI, phi2)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) orimat(3,3), phi1, PHI, phi2
 real(pReal) sqhkl, squvw, sqhk, val

 sqhkl=sqrt(orimat(1,3)*orimat(1,3)+orimat(2,3)*orimat(2,3)+orimat(3,3)*orimat(3,3))
 squvw=sqrt(orimat(1,1)*orimat(1,1)+orimat(2,1)*orimat(2,1)+orimat(3,1)*orimat(3,1))
 sqhk=sqrt(orimat(1,3)*orimat(1,3)+orimat(2,3)*orimat(2,3))
! calculate PHI
 val=orimat(3,3)/sqhkl
 
 if(val.GT.1.0_pReal) val=1.0_pReal
 if(val.LT.-1.0_pReal) val=-1.0_pReal
     
 PHI=acos(val)

 if(PHI.LT.1.0e-30_pReal) then
! calculate phi2
     phi2=0.0
! calculate phi1
     val=orimat(1,1)/squvw
 
     if(val.GT.1.0_pReal) val=1.0_pReal
     if(val.LT.-1.0_pReal) val=-1.0_pReal
     
     if(orimat(2,1).LE.0.0) then
    phi1=acos(val)
     else
    phi1=2.0_pReal*pi-acos(val)
     end if
 else
! calculate phi2
     val=orimat(2,3)/sqhk
 
     if(val.GT.1.0_pReal) val=1.0_pReal
     if(val.LT.-1.0_pReal) val=-1.0_pReal
     
     if(orimat(1,3).GE.0.0) then
    phi2=acos(val)
     else
    phi2=2.0_pReal*pi-acos(val)
     end if
! calculate phi1
     val=-orimat(3,2)/sin(PHI)
 
     if(val.GT.1.0_pReal) val=1.0_pReal
     if(val.LT.-1.0_pReal) val=-1.0_pReal
     
     if(orimat(3,1).GE.0.0) then
    phi1=acos(val)
     else
    phi1=2.0_pReal*pi-acos(val)
     end if
 end if
! convert angles to degrees
 phi1=phi1*inDeg
 PHI=PHI*inDeg
 phi2=phi2*inDeg
 return
 
 END SUBROUTINE


!#####################################################
! bestimmt Drehmatrix DREH3 fuer Drehung um Omega um Achse (u,v,w)
 FUNCTION math_RodrigtoR(Omega,U,V,W)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) omega, u, v, w, math_RodrigtoR(3,3)
 real(pReal) betrag, s, c, u2, v2, w2

 BETRAG=SQRT(U**2+V**2+W**2)
 S=SIN(OMEGA)
 C=COS(OMEGA)
 U2=U/BETRAG
 V2=V/BETRAG
 W2=W/BETRAG
 math_RodrigtoR(1,1)=(1-U2**2)*C+U2**2
 math_RodrigtoR(1,2)=U2*V2*(1-C)+W2*S
 math_RodrigtoR(1,3)=U2*W2*(1-C)-V2*S
 math_RodrigtoR(2,1)=U2*V2*(1-C)-W2*S
 math_RodrigtoR(2,2)=(1-V2**2)*C+V2**2
 math_RodrigtoR(2,3)=V2*W2*(1-C)+U2*S
 math_RodrigtoR(3,1)=U2*W2*(1-C)+V2*S
 math_RodrigtoR(3,2)=V2*W2*(1-C)-U2*S
 math_RodrigtoR(3,3)=(1-W2**2)*C+W2**2
 return

 END FUNCTION

! Best. Drehmatrix ROTA fuer Euler-Winkel  

 FUNCTION math_EulertoR (P1,P,P2)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) p1, p, p2, math_EulertoR(3,3)
 real(pReal)  xp1, xp, xp2, c1, c, c2, s1, s, s2

 XP1=P1*inRad
 XP=P*inRad
 XP2=P2*inRad
 C1=COS(XP1)
 C=COS(XP)
 C2=COS(XP2)
 S1=SIN(XP1)
 S=SIN(XP)
 S2=SIN(XP2)
 math_EulertoR(1,1)=C1*C2-S1*S2*C
 math_EulertoR(1,2)=S1*C2+C1*S2*C
 math_EulertoR(1,3)=S2*S
 math_EulertoR(2,1)=-C1*S2-S1*C2*C
 math_EulertoR(2,2)=-S1*S2+C1*C2*C
 math_EulertoR(2,3)=C2*S
 math_EulertoR(3,1)=S1*S
 math_EulertoR(3,2)=-C1*S
 math_EulertoR(3,3)=C
 return
 
 END FUNCTION


!**************************************************************************
! BERECHNUNG VON ORIENTIERUNGSBEZIEHUNGEN ZWISCHEN
! ZWEI VORGEGEBENEN ORIENTIERUNGEN

 function math_disorient(P1,P,P2)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) D1(3,3),D2(3,3),P1(2),P(2),P2(2),D1T(3,3),DR(3,3)
 real(pReal) math_disorient, spur, sp, omega, alpha
 integer(pInt) i

! ERSTELLEN DER BEIDEN DMATRIZEN

 d1 = math_EulertoR(p1(1),P(1),p2(1))
 d2 = math_EulertoR(p1(2),P(2),p2(2))
!****************************************************
! BESTIMMUNG DER INVERSEN MATRIX ZUR ORIENTIERUNG 1:DM
!****************************************************
 d1T=transpose(d1)
!***********************************************************
! MATRIZENMULTIPLIKATION DER MATRIZEN D2 UND DM=DR(I,J)
!***********************************************************
 dr=matmul(d2,d1T)
!*******************************
! BESTIMMUNG DES ROTATIONSWINKELS
!*******************************
 SPUR=DR(1,1)+DR(2,2)+DR(3,3)
 SP=(SPUR-1._pReal)*0.4999999_pReal
 OMEGA=PI*0.5_pReal-ASIN(SP)
! Winkel in Grad umrechnen
 ALPHA=OMEGA*inDeg
 math_disorient=abs(alpha)
 return

 END FUNCTION


!****************************************************************
 subroutine math_pDecomposition(FE,U,R,ISING)
!-----FE=RU 
!-----INVERT is the subroutine applied by Marc
!****************************************************************
 use prec, only: pReal, pInt
 implicit none

 integer(pInt)  ISING
 real(pReal) FE(3,3),R(3,3),U(3,3),CE(3,3),EW1,EW2,EW3,EB1(3,3),EB2(3,3),EB3(3,3),UI(3,3),det
 ising=0
 ce=matmul(transpose(fe),fe)
 CALL math_spectral1(CE,EW1,EW2,EW3,EB1,EB2,EB3)
 U=DSQRT(EW1)*EB1+DSQRT(EW2)*EB2+DSQRT(EW3)*EB3
 UI=U
 call invert(UI,3,0,0,det,3)
 if (det.EQ.0) then
     ising=1
     return
 endif
 R=matmul(fe,ui)
 return 
 
 END SUBROUTINE


!**********************************************************************
 subroutine math_spectral1(M,EW1,EW2,EW3,EB1,EB2,EB3)
!**** EIGENWERTE UND EIGENWERTBASIS DER SYMMETRISCHEN 3X3 MATRIX M

 use prec, only: pReal, pInt
 implicit none

 real(pReal) M(3,3),EB1(3,3),EB2(3,3),EB3(3,3),EW1,EW2,EW3
 real(pReal) HI1M,HI2M,HI3M,TOL,R,S,T,P,Q,RHO,PHI,Y1,Y2,Y3,D1,D2,D3
 real(pReal) C1,C2,C3,M1(3,3),M2(3,3),M3(3,3),I3(3,3),arg
 TOL=1.e-14_pReal
 CALL math_hi(M,HI1M,HI2M,HI3M)
 R=-HI1M
 S= HI2M
 T=-HI3M
 P=S-R**2.0_pReal/3.0_pReal
 Q=2.0_pReal/27.0_pReal*R**3.0_pReal-R*S/3.0_pReal+T
 EB1=0.0_pReal
 EB2=0.0_pReal
 EB3=0.0_pReal
 IF((ABS(P).LT.TOL).AND.(ABS(Q).LT.TOL))THEN
!   DREI GLEICHE EIGENWERTE
   EW1=HI1M/3.0_pReal
   EW2=EW1
   EW3=EW1
!   this is not really correct, but this way U is calculated
!   correctly in PDECOMPOSITION (correct is EB?=I)
   EB1(1,1)=1.0_pReal
   EB2(2,2)=1.0_pReal
   EB3(3,3)=1.0_pReal
 ELSE
   RHO=SQRT(-3.0_pReal*P**3.0_pReal)/9.0_pReal
   arg=-Q/RHO/2.0_pReal
   if(arg.GT.1) arg=1
   if(arg.LT.-1) arg=-1
   PHI=ACOS(arg)
   Y1=2*RHO**(1.0_pReal/3.0_pReal)*COS(PHI/3.0_pReal)
   Y2=2*RHO**(1.0_pReal/3.0_pReal)*COS(PHI/3.0_pReal+2.0_pReal/3.0_pReal*PI)
   Y3=2*RHO**(1.0_pReal/3.0_pReal)*COS(PHI/3.0_pReal+4.0_pReal/3.0_pReal*PI)
   EW1=Y1-R/3.0_pReal
   EW2=Y2-R/3.0_pReal
   EW3=Y3-R/3.0_pReal
   C1=ABS(EW1-EW2)
   C2=ABS(EW2-EW3) 
   C3=ABS(EW3-EW1)

   IF(C1.LT.TOL) THEN
!  EW1 is equal to EW2
  D3=1.0_pReal/(EW3-EW1)/(EW3-EW2)
  M1=M-EW1*math_I3
  M2=M-EW2*math_I3
  EB3=MATMUL(M1,M2)*D3
  EB1=math_I3-EB3
!  both EB2 and EW2 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW2=0.0_pReal
   ELSE IF(C2.LT.TOL) THEN
!  EW2 is equal to EW3
  D1=1.0_pReal/(EW1-EW2)/(EW1-EW3)
  M2=M-math_I3*EW2
  M3=M-math_I3*EW3
  EB1=MATMUL(M2,M3)*D1
  EB2=math_I3-EB1
!  both EB3 and EW3 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW3=0.0_pReal
   ELSE IF(C3.LT.TOL) THEN
!  EW1 is equal to EW3
  D2=1.0_pReal/(EW2-EW1)/(EW2-EW3) 
  M1=M-math_I3*EW1
  M3=M-math_I3*EW3
  EB2=MATMUL(M1,M3)*D2
  EB1=math_I3-EB2
!  both EB3 and EW3 are set to zero so that they do not
!  contribute to U in PDECOMPOSITION
  EW3=0.0_pReal
   ELSE
!  all three eigenvectors are different
  D1=1.0_pReal/(EW1-EW2)/(EW1-EW3)
  D2=1.0_pReal/(EW2-EW1)/(EW2-EW3) 
  D3=1.0_pReal/(EW3-EW1)/(EW3-EW2)
  M1=M-EW1*math_I3
  M2=M-EW2*math_I3
  M3=M-EW3*math_I3
  EB1=MATMUL(M2,M3)*D1
  EB2=MATMUL(M1,M3)*D2
  EB3=MATMUL(M1,M2)*D3
   END IF
 END IF
 RETURN
 END SUBROUTINE


!**********************************************************************
!**** EINHEITSMATRIX MIT dim DIAGONALELEMENTEN 

 FUNCTION math_identity2nd(dimen)  

 use prec, only: pReal, pInt
 implicit none

 integer(pInt)  i,dimen
 real(pReal) math_identity2nd(dimen,dimen)

 math_identity2nd = 0.0_pReal 
 forall (i=1:dimen) math_identity2nd(i,i) = 1.0_pReal 
 return

 END FUNCTION

!**********************************************************************
!**** EINHEITSTENSOR 4th MIT dim "DIAGONAL"ELEMENTEN 

 FUNCTION math_identity4th(dimen)  

 use prec, only: pReal, pInt
 implicit none

 integer(pInt)  i,j,k,l,dimen
 real(pReal) math_identity4th(dimen,dimen,dimen,dimen)

 forall (i=1:dimen,j=1:dimen,k=1:dimen,l=1:dimen) math_identity4th(i,j,k,l) = &
 	0.5_pReal*(math_I3(i,k)*math_I3(j,k)+math_I3(i,l)*math_I3(j,k)) 
 return

 END FUNCTION


!********************************************************************** 
!**** HAUPTINVARIANTEN HI1M, HI2M, HI3M DER 3X3 MATRIX M

 SUBROUTINE math_hi(M,HI1M,HI2M,HI3M)
 use prec, only: pReal, pInt
 implicit none

 real(pReal) M(3,3),HI1M,HI2M,HI3M 

 HI1M=M(1,1)+M(2,2)+M(3,3)
 HI2M=(M(1,1)+M(2,2)+M(3,3))**2/2.0_pReal-(M(1,1)**2+M(2,2)**2+M(3,3)**2)/2.0_pReal-M(1,2)*M(2,1)-M(1,3)*M(3,1)-M(2,3)*M(3,2)
 HI3M=math_det3x3(M)
! QUESTION: is 3rd equiv det(M) ?? if yes, use function math_det !agreed on YES
 return  

 END SUBROUTINE


 SUBROUTINE get_seed(seed)
!
!*******************************************************************************
!
!! GET_SEED returns a seed for the random number generator.
!
!
!  Discussion:
!
! The seed depends on the current time, and ought to be (slightly)
! different every millisecond. Once the seed is obtained, a random
! number generator should be called a few times to further process
! the seed.
!
!  Modified:
!
! 27 June 2000
!
!  Author:
!
! John Burkardt
!
!  Parameters:
!
! Output, integer SEED, a pseudorandom seed value.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, only: pReal, pInt
 implicit none

 integer(pInt) seed
 real(pReal)  temp
 character ( len = 10 ) time
 character ( len = 8 ) today
 integer(pInt) values(8)
 character ( len = 5 ) zone

 call date_and_time ( today, time, zone, values )

 temp = 0.0D+00

 temp = temp + dble ( values(2) - 1 ) / 11.0D+00
 temp = temp + dble ( values(3) - 1 ) / 30.0D+00
 temp = temp + dble ( values(5) ) / 23.0D+00
 temp = temp + dble ( values(6) ) / 59.0D+00
 temp = temp + dble ( values(7) ) / 59.0D+00
 temp = temp + dble ( values(8) ) / 999.0D+00
 temp = temp / 6.0D+00

 if ( temp <= 0.0D+00 ) then
   temp = 1.0D+00 / 3.0D+00
 else if ( 1.0D+00 <= temp ) then
   temp = 2.0D+00 / 3.0D+00
 end if

 seed = int ( dble ( huge ( 1 ) ) * temp , pInt)
!
!  Never use a seed of 0 or maximum integer.
!
 if ( seed == 0 ) then
   seed = 1
 end if

 if ( seed == huge ( 1 ) ) then
   seed = seed - 1
 end if

 return

 END SUBROUTINE


 subroutine halton ( ndim, r )
!
!*******************************************************************************
!
!! HALTON computes the next element in the Halton sequence.
!
!
!  Modified:
!
! 09 March 2003
!
!  Author:
!
! John Burkardt
!
!  Parameters:
!
! Input, integer NDIM, the dimension of the element.
!
! Output, real R(NDIM), the next element of the current Halton
! sequence.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, ONLY: pReal, pInt
 implicit none

 integer(pInt) ndim

 integer(pInt) base(ndim)
 real(pReal) r(ndim)
 integer(pInt) seed
 integer(pInt) value(1)

 call halton_memory ( 'GET', 'SEED', 1, value )
 seed = value(1)

 call halton_memory ( 'GET', 'BASE', ndim, base )

 call i_to_halton ( seed, base, ndim, r )

 value(1) = 1
 call halton_memory ( 'INC', 'SEED', 1, value )

 return

 END SUBROUTINE


 subroutine halton_memory ( action, name, ndim, value )
!
!*******************************************************************************
!
!! HALTON_MEMORY sets or returns quantities associated with the Halton sequence.
!
!
!  Modified:
!
! 09 March 2003
!
!  Author:
!
! John Burkardt
!
!  Parameters:
!
! Input, character ( len = * ) ACTION, the desired action.
! 'GET' means get the value of a particular quantity.
! 'SET' means set the value of a particular quantity.
! 'INC' means increment the value of a particular quantity.
!  (Only the SEED can be incremented.)
!
! Input, character ( len = * ) NAME, the name of the quantity.
! 'BASE' means the Halton base or bases.
! 'NDIM' means the spatial dimension.
! 'SEED' means the current Halton seed.
!
! Input/output, integer NDIM, the dimension of the quantity.
! If ACTION is 'SET' and NAME is 'BASE', then NDIM is input, and
! is the number of entries in VALUE to be put into BASE.
!
! Input/output, integer VALUE(NDIM), contains a value.
! If ACTION is 'SET', then on input, VALUE contains values to be assigned
! to the internal variable.
! If ACTION is 'GET', then on output, VALUE contains the values of
! the specified internal variable.
! If ACTION is 'INC', then on input, VALUE contains the increment to
! be added to the specified internal variable.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, only: pReal, pInt
 implicit none

 character ( len = * ) action
 integer(pInt), allocatable, save :: base(:)
 logical, save :: first_call = .true.
 integer(pInt) i
 character ( len = * ) name
 integer(pInt) ndim
 integer(pInt), save :: ndim_save = 0
 integer(pInt) prime
 integer(pInt), save :: seed = 1
 integer(pInt) value(*)

 if ( first_call ) then
   ndim_save = 1
   allocate ( base(ndim_save) )
   base(1) = 2
   first_call = .false.
 end if
!
!  Set
!
 if ( action(1:1) == 'S' .or. action(1:1) == 's' ) then

   if ( name(1:1) == 'B' .or. name(1:1) == 'b' ) then

     if ( ndim_save /= ndim ) then
  deallocate ( base )
  ndim_save = ndim
  allocate ( base(ndim_save) )
     end if

     base(1:ndim) = value(1:ndim)

   else if ( name(1:1) == 'N' .or. name(1:1) == 'n' ) then

     if ( ndim_save /= value(1) ) then
  deallocate ( base )
  ndim_save = value(1)
  allocate ( base(ndim_save) )
  do i = 1, ndim_save
    base(i) = prime ( i )
  end do
     else
  ndim_save = value(1)
     end if

   else if ( name(1:1) == 'S' .or. name(1:1) == 's' ) then

     seed = value(1)

 end if
!
!  Get
!
 else if ( action(1:1) == 'G' .or. action(1:1) == 'g' ) then

   if ( name(1:1) == 'B' .or. name(1:1) == 'b' ) then

     if ( ndim /= ndim_save ) then
  deallocate ( base )
  ndim_save = ndim
  allocate ( base(ndim_save) )
  do i = 1, ndim_save
    base(i) = prime(i)
  end do
     end if

     value(1:ndim_save) = base(1:ndim_save)

   else if ( name(1:1) == 'N' .or. name(1:1) == 'n' ) then

     value(1) = ndim_save

   else if ( name(1:1) == 'S' .or. name(1:1) == 's' ) then

     value(1) = seed

   end if
!
!  Increment
!
 else if ( action(1:1) == 'I' .or. action(1:1) == 'i' ) then

   if ( name(1:1) == 'S' .or. name(1:1) == 's' ) then
     seed = seed + value(1)
   end if

 end if

 return

 END SUBROUTINE



 subroutine halton_ndim_set ( ndim )
!
!*******************************************************************************
!
!! HALTON_NDIM_SET sets the dimension for a Halton sequence.
!
!
!  Modified:
!
! 26 February 2001
!
!  Author:
!
! John Burkardt
!
!  Parameters:
!
! Input, integer NDIM, the dimension of the Halton vectors.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, only: pReal, pInt
 implicit none

 integer(pInt) ndim
 integer(pInt) value(1)

 value(1) = ndim
 call halton_memory ( 'SET', 'NDIM', 1, value )

 return

 END SUBROUTINE


 subroutine halton_seed_set ( seed )
!
!*******************************************************************************
!
!! HALTON_SEED_SET sets the "seed" for the Halton sequence.
!
!
!  Discussion:
!
! Calling HALTON repeatedly returns the elements of the
! Halton sequence in order, starting with element number 1.
! An internal counter, called SEED, keeps track of the next element
! to return. Each time the routine is called, the SEED-th element
! is computed, and then SEED is incremented by 1.
!
! To restart the Halton sequence, it is only necessary to reset
! SEED to 1. It might also be desirable to reset SEED to some other value.
! This routine allows the user to specify any value of SEED.
!
! The default value of SEED is 1, which restarts the Halton sequence.
!
!  Modified:
!
! 26 February 2001
!
!  Author:
!
! John Burkardt
!
!  Parameters:
!
! Input, integer SEED, the seed for the Halton sequence.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, only: pReal, pInt
 implicit none

 integer(pInt), parameter :: ndim = 1

 integer(pInt) seed
 integer(pInt) value(ndim)

 value(1) = seed
 call halton_memory ( 'SET', 'SEED', ndim, value )

 return

 END SUBROUTINE


 subroutine i_to_halton ( seed, base, ndim, r )
!
!*******************************************************************************
!
!! I_TO_HALTON computes an element of a Halton sequence.
!
!
!  Reference:
!
! J H Halton,
! On the efficiency of certain quasi-random sequences of points
! in evaluating multi-dimensional integrals,
! Numerische Mathematik,
! Volume 2, pages 84-90, 1960.
!
!  Modified:
!
! 26 February 2001
!
!  Author:
!
! John Burkardt
!
!  Parameters:
!
! Input, integer SEED, the index of the desired element.
! Only the absolute value of SEED is considered. SEED = 0 is allowed,
! and returns R = 0.
!
! Input, integer BASE(NDIM), the Halton bases, which should be
! distinct prime numbers.  This routine only checks that each base
! is greater than 1.
!
! Input, integer NDIM, the dimension of the sequence.
!
! Output, real R(NDIM), the SEED-th element of the Halton sequence
! for the given bases.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, ONLY: pReal, pInt
 implicit none

 integer(pInt) ndim

 integer(pInt) base(ndim)
 real(pReal) base_inv(ndim)
 integer(pInt) digit(ndim)
 integer(pInt) i
 real(pReal) r(ndim)
 integer(pInt) seed
 integer(pInt) seed2(ndim)

 seed2(1:ndim) = abs ( seed )

 r(1:ndim) = 0.0_pReal

 if ( any ( base(1:ndim) <= 1 ) ) then
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'I_TO_HALTON - Fatal error!'
   write ( *, '(a)' ) ' An input base BASE is <= 1!'
   do i = 1, ndim
     write ( *, '(i6,i6)' ) i, base(i)
   end do
   call flush(6)
   stop
 end if

 base_inv(1:ndim) = 1.0_pReal / real ( base(1:ndim), pReal )

 do while ( any ( seed2(1:ndim) /= 0 ) )
   digit(1:ndim) = mod ( seed2(1:ndim), base(1:ndim) )
   r(1:ndim) = r(1:ndim) + real ( digit(1:ndim), pReal ) * base_inv(1:ndim)
   base_inv(1:ndim) = base_inv(1:ndim) / real ( base(1:ndim), pReal )
   seed2(1:ndim) = seed2(1:ndim) / base(1:ndim)
 end do

 return

 END SUBROUTINE


 function prime ( n )
!
!*******************************************************************************
!
!! PRIME returns any of the first PRIME_MAX prime numbers.
!
!
!  Note:
!
! PRIME_MAX is 1500, and the largest prime stored is 12553.
!
!  Modified:
!
! 21 June 2002
!
!  Author:
!
! John Burkardt
!
!  Reference:
!
! Milton Abramowitz and Irene Stegun,
! Handbook of Mathematical Functions,
! US Department of Commerce, 1964, pages 870-873.
!
! Daniel Zwillinger,
! CRC Standard Mathematical Tables and Formulae,
! 30th Edition,
! CRC Press, 1996, pages 95-98.
!
!  Parameters:
!
! Input, integer N, the index of the desired prime number.
! N = -1 returns PRIME_MAX, the index of the largest prime available.
! N = 0 is legal, returning PRIME = 1.
! It should generally be true that 0 <= N <= PRIME_MAX.
!
! Output, integer PRIME, the N-th prime.  If N is out of range, PRIME
! is returned as 0.
!
!  Modified:
!
! 29 April 2005
!
!  Author:
!
! Franz Roters
!
 use prec, only: pReal, pInt
 implicit none

 integer(pInt), parameter :: prime_max = 1500

 integer(pInt), save :: icall = 0
 integer(pInt) n
 integer(pInt), save, dimension ( prime_max ) :: npvec
 integer(pInt) prime

 if ( icall == 0 ) then

 icall = 1

 npvec(1:100) = (/&
        2,    3,    5,    7,   11,   13,   17,   19,   23,   29, &
       31,   37,   41,   43,   47,   53,   59,   61,   67,   71, &
       73,   79,   83,   89,   97,  101,  103,  107,  109,  113, &
      127,  131,  137,  139,  149,  151,  157,  163,  167,  173, &
      179,  181,  191,  193,  197,  199,  211,  223,  227,  229, &
      233,  239,  241,  251,  257,  263,  269,  271,  277,  281, &
      283,  293,  307,  311,  313,  317,  331,  337,  347,  349, &
      353,  359,  367,  373,  379,  383,  389,  397,  401,  409, &
      419,  421,  431,  433,  439,  443,  449,  457,  461,  463, &
      467,  479,  487,  491,  499,  503,  509,  521,  523,  541 /)

 npvec(101:200) = (/ &
       547,  557,  563,  569,  571,  577,  587,  593,  599,  601, &
       607,  613,  617,  619,  631,  641,  643,  647,  653,  659, &
       661,  673,  677,  683,  691,  701,  709,  719,  727,  733, &
       739,  743,  751,  757,  761,  769,  773,  787,  797,  809, &
       811,  821,  823,  827,  829,  839,  853,  857,  859,  863, &
       877,  881,  883,  887,  907,  911,  919,  929,  937,  941, &
       947,  953,  967,  971,  977,  983,  991,  997, 1009, 1013, &
      1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
      1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
      1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223 /)

 npvec(201:300) = (/ &
      1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
      1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
      1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
      1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
      1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
      1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
      1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
      1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
      91823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
      1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987 /)

 npvec(301:400) = (/ &
      1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053, &
      2063, 2069, 2081, 2083, 2087, 2089, 2099, 2111, 2113, 2129, &
      2131, 2137, 2141, 2143, 2153, 2161, 2179, 2203, 2207, 2213, &
      2221, 2237, 2239, 2243, 2251, 2267, 2269, 2273, 2281, 2287, &
      2293, 2297, 2309, 2311, 2333, 2339, 2341, 2347, 2351, 2357, &
      2371, 2377, 2381, 2383, 2389, 2393, 2399, 2411, 2417, 2423, &
      2437, 2441, 2447, 2459, 2467, 2473, 2477, 2503, 2521, 2531, &
      2539, 2543, 2549, 2551, 2557, 2579, 2591, 2593, 2609, 2617, &
      2621, 2633, 2647, 2657, 2659, 2663, 2671, 2677, 2683, 2687, &
      2689, 2693, 2699, 2707, 2711, 2713, 2719, 2729, 2731, 2741 /)

 npvec(401:500) = (/ &
      2749, 2753, 2767, 2777, 2789, 2791, 2797, 2801, 2803, 2819, &
      2833, 2837, 2843, 2851, 2857, 2861, 2879, 2887, 2897, 2903, &
      2909, 2917, 2927, 2939, 2953, 2957, 2963, 2969, 2971, 2999, &
      3001, 3011, 3019, 3023, 3037, 3041, 3049, 3061, 3067, 3079, &
      3083, 3089, 3109, 3119, 3121, 3137, 3163, 3167, 3169, 3181, &
      3187, 3191, 3203, 3209, 3217, 3221, 3229, 3251, 3253, 3257, &
      3259, 3271, 3299, 3301, 3307, 3313, 3319, 3323, 3329, 3331, &
      3343, 3347, 3359, 3361, 3371, 3373, 3389, 3391, 3407, 3413, &
      3433, 3449, 3457, 3461, 3463, 3467, 3469, 3491, 3499, 3511, &
      3517, 3527, 3529, 3533, 3539, 3541, 3547, 3557, 3559, 3571 /)

 npvec(501:600) = (/ &
      3581, 3583, 3593, 3607, 3613, 3617, 3623, 3631, 3637, 3643, &
      3659, 3671, 3673, 3677, 3691, 3697, 3701, 3709, 3719, 3727, &
      3733, 3739, 3761, 3767, 3769, 3779, 3793, 3797, 3803, 3821, &
      3823, 3833, 3847, 3851, 3853, 3863, 3877, 3881, 3889, 3907, &
      3911, 3917, 3919, 3923, 3929, 3931, 3943, 3947, 3967, 3989, &
      4001, 4003, 4007, 4013, 4019, 4021, 4027, 4049, 4051, 4057, &
      4073, 4079, 4091, 4093, 4099, 4111, 4127, 4129, 4133, 4139, &
      4153, 4157, 4159, 4177, 4201, 4211, 4217, 4219, 4229, 4231, &
      4241, 4243, 4253, 4259, 4261, 4271, 4273, 4283, 4289, 4297, &
      4327, 4337, 4339, 4349, 4357, 4363, 4373, 4391, 4397, 4409 /)

 npvec(601:700) = (/ &
      4421, 4423, 4441, 4447, 4451, 4457, 4463, 4481, 4483, 4493, &
      4507, 4513, 4517, 4519, 4523, 4547, 4549, 4561, 4567, 4583, &
      4591, 4597, 4603, 4621, 4637, 4639, 4643, 4649, 4651, 4657, &
      4663, 4673, 4679, 4691, 4703, 4721, 4723, 4729, 4733, 4751, &
      4759, 4783, 4787, 4789, 4793, 4799, 4801, 4813, 4817, 4831, &
      4861, 4871, 4877, 4889, 4903, 4909, 4919, 4931, 4933, 4937, &
      4943, 4951, 4957, 4967, 4969, 4973, 4987, 4993, 4999, 5003, &
      5009, 5011, 5021, 5023, 5039, 5051, 5059, 5077, 5081, 5087, &
      5099, 5101, 5107, 5113, 5119, 5147, 5153, 5167, 5171, 5179, &
      5189, 5197, 5209, 5227, 5231, 5233, 5237, 5261, 5273, 5279 /)

 npvec(701:800) = (/ &
      5281, 5297, 5303, 5309, 5323, 5333, 5347, 5351, 5381, 5387, &
      5393, 5399, 5407, 5413, 5417, 5419, 5431, 5437, 5441, 5443, &
      5449, 5471, 5477, 5479, 5483, 5501, 5503, 5507, 5519, 5521, &
      5527, 5531, 5557, 5563, 5569, 5573, 5581, 5591, 5623, 5639, &
      5641, 5647, 5651, 5653, 5657, 5659, 5669, 5683, 5689, 5693, &
      5701, 5711, 5717, 5737, 5741, 5743, 5749, 5779, 5783, 5791, &
      5801, 5807, 5813, 5821, 5827, 5839, 5843, 5849, 5851, 5857, &
      5861, 5867, 5869, 5879, 5881, 5897, 5903, 5923, 5927, 5939, &
      5953, 5981, 5987, 6007, 6011, 6029, 6037, 6043, 6047, 6053, &
      6067, 6073, 6079, 6089, 6091, 6101, 6113, 6121, 6131, 6133 /)

 npvec(801:900) = (/ &
      6143, 6151, 6163, 6173, 6197, 6199, 6203, 6211, 6217, 6221, &
      6229, 6247, 6257, 6263, 6269, 6271, 6277, 6287, 6299, 6301, &
      6311, 6317, 6323, 6329, 6337, 6343, 6353, 6359, 6361, 6367, &
      6373, 6379, 6389, 6397, 6421, 6427, 6449, 6451, 6469, 6473, &
      6481, 6491, 6521, 6529, 6547, 6551, 6553, 6563, 6569, 6571, &
      6577, 6581, 6599, 6607, 6619, 6637, 6653, 6659, 6661, 6673, &
      6679, 6689, 6691, 6701, 6703, 6709, 6719, 6733, 6737, 6761, &
      6763, 6779, 6781, 6791, 6793, 6803, 6823, 6827, 6829, 6833, &
      6841, 6857, 6863, 6869, 6871, 6883, 6899, 6907, 6911, 6917, &
      6947, 6949, 6959, 6961, 6967, 6971, 6977, 6983, 6991, 6997 /)

 npvec(901:1000) = (/ &
      7001, 7013, 7019, 7027, 7039, 7043, 7057, 7069, 7079, 7103, &
      7109, 7121, 7127, 7129, 7151, 7159, 7177, 7187, 7193, 7207, &
      7211, 7213, 7219, 7229, 7237, 7243, 7247, 7253, 7283, 7297, &
      7307, 7309, 7321, 7331, 7333, 7349, 7351, 7369, 7393, 7411, &
      7417, 7433, 7451, 7457, 7459, 7477, 7481, 7487, 7489, 7499, &
      7507, 7517, 7523, 7529, 7537, 7541, 7547, 7549, 7559, 7561, &
      7573, 7577, 7583, 7589, 7591, 7603, 7607, 7621, 7639, 7643, &
      7649, 7669, 7673, 7681, 7687, 7691, 7699, 7703, 7717, 7723, &
      7727, 7741, 7753, 7757, 7759, 7789, 7793, 7817, 7823, 7829, &
      7841, 7853, 7867, 7873, 7877, 7879, 7883, 7901, 7907, 7919 /)

 npvec(1001:1100) = (/ &
      7927, 7933, 7937, 7949, 7951, 7963, 7993, 8009, 8011, 8017, &
      8039, 8053, 8059, 8069, 8081, 8087, 8089, 8093, 8101, 8111, &
      8117, 8123, 8147, 8161, 8167, 8171, 8179, 8191, 8209, 8219, &
      8221, 8231, 8233, 8237, 8243, 8263, 8269, 8273, 8287, 8291, &
      8293, 8297, 8311, 8317, 8329, 8353, 8363, 8369, 8377, 8387, &
      8389, 8419, 8423, 8429, 8431, 8443, 8447, 8461, 8467, 8501, &
      8513, 8521, 8527, 8537, 8539, 8543, 8563, 8573, 8581, 8597, &
      8599, 8609, 8623, 8627, 8629, 8641, 8647, 8663, 8669, 8677, &
      8681, 8689, 8693, 8699, 8707, 8713, 8719, 8731, 8737, 8741, &
      8747, 8753, 8761, 8779, 8783, 8803, 8807, 8819, 8821, 8831 /)

 npvec(1101:1200) = (/ &
      8837, 8839, 8849, 8861, 8863, 8867, 8887, 8893, 8923, 8929, &
      8933, 8941, 8951, 8963, 8969, 8971, 8999, 9001, 9007, 9011, &
      9013, 9029, 9041, 9043, 9049, 9059, 9067, 9091, 9103, 9109, &
      9127, 9133, 9137, 9151, 9157, 9161, 9173, 9181, 9187, 9199, &
      9203, 9209, 9221, 9227, 9239, 9241, 9257, 9277, 9281, 9283, &
      9293, 9311, 9319, 9323, 9337, 9341, 9343, 9349, 9371, 9377, &
      9391, 9397, 9403, 9413, 9419, 9421, 9431, 9433, 9437, 9439, &
      9461, 9463, 9467, 9473, 9479, 9491, 9497, 9511, 9521, 9533, &
      9539, 9547, 9551, 9587, 9601, 9613, 9619, 9623, 9629, 9631, &
      9643, 9649, 9661, 9677, 9679, 9689, 9697, 9719, 9721, 9733 /)

 npvec(1201:1300) = (/ &
       9739, 9743, 9749, 9767, 9769, 9781, 9787, 9791, 9803, 9811, &
       9817, 9829, 9833, 9839, 9851, 9857, 9859, 9871, 9883, 9887, &
       9901, 9907, 9923, 9929, 9931, 9941, 9949, 9967, 9973,10007, &
      10009,10037,10039,10061,10067,10069,10079,10091,10093,10099, &
      10103,10111,10133,10139,10141,10151,10159,10163,10169,10177, &
      10181,10193,10211,10223,10243,10247,10253,10259,10267,10271, &
      10273,10289,10301,10303,10313,10321,10331,10333,10337,10343, &
      10357,10369,10391,10399,10427,10429,10433,10453,10457,10459, &
      10463,10477,10487,10499,10501,10513,10529,10531,10559,10567, &
      10589,10597,10601,10607,10613,10627,10631,10639,10651,10657 /)

 npvec(1301:1400) = (/ &
      10663,10667,10687,10691,10709,10711,10723,10729,10733,10739, &
      10753,10771,10781,10789,10799,10831,10837,10847,10853,10859, &
      10861,10867,10883,10889,10891,10903,10909,19037,10939,10949, &
      10957,10973,10979,10987,10993,11003,11027,11047,11057,11059, &
      11069,11071,11083,11087,11093,11113,11117,11119,11131,11149, &
      11159,11161,11171,11173,11177,11197,11213,11239,11243,11251, &
      11257,11261,11273,11279,11287,11299,11311,11317,11321,11329, &
      11351,11353,11369,11383,11393,11399,11411,11423,11437,11443, &
      11447,11467,11471,11483,11489,11491,11497,11503,11519,11527, &
      11549,11551,11579,11587,11593,11597,11617,11621,11633,11657 /)

 npvec(1401:1500) = (/ &
      11677,11681,11689,11699,11701,11717,11719,11731,11743,11777, &
      11779,11783,11789,11801,11807,11813,11821,11827,11831,11833, &
      11839,11863,11867,11887,11897,11903,11909,11923,11927,11933, &
      11939,11941,11953,11959,11969,11971,11981,11987,12007,12011, &
      12037,12041,12043,12049,12071,12073,12097,12101,12107,12109, &
      12113,12119,12143,12149,12157,12161,12163,12197,12203,12211, &
      12227,12239,12241,12251,12253,12263,12269,12277,12281,12289, &
      12301,12323,12329,12343,12347,12373,12377,12379,12391,12401, &
      12409,12413,12421,12433,12437,12451,12457,12473,12479,12487, &
      12491,12497,12503,12511,12517,12527,12539,12541,12547,12553 /)

 end if

 if ( n == -1 ) then
   prime = prime_max
 else if ( n == 0 ) then
   prime = 1
 else if ( n <= prime_max ) then
   prime = npvec(n)
 else
   prime = 0
   write ( *, '(a)' ) ' '
   write ( *, '(a)' ) 'PRIME - Fatal error!'
   write ( *, '(a,i6)' ) '  Illegal prime index N = ', n
   write ( *, '(a,i6)' ) '  N must be between 0 and PRIME_MAX =',prime_max
   call flush(6)
   stop
 end if

 return
 
 END FUNCTION


!********************************************************************
! This routine generates a random orientation
!********************************************************************
 subroutine math_random_ori (phi1, PHI, phi2, scatter)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) phi1, PHI, phi2, scatter, x, y, z

 call random_number(x)
 call random_number(y)
 call random_number(z)
 phi1=x*360.0_pReal
 PHI=acos(y)*inDeg
 phi2=z*360.0_pReal
 scatter=0.0_pReal
 return

 END SUBROUTINE


 subroutine math_halton_ori (phi1, PHI, phi2, scatter)
!********************************************************************
! This routine generates a random orientation using Halton series
!********************************************************************
 use prec, only: pReal, pInt
 implicit none

 real(pReal) phi1, PHI, phi2, scatter, r(3)

 call halton(3,r)
 phi1=r(1)*360.0_pReal
 PHI=acos(r(2))*inDeg
 phi2=r(3)*360.0_pReal
 scatter=0.0_pReal
 return
 
 END SUBROUTINE


!********************************************************************
! This routine applies gaussian scatter to the texture components
!********************************************************************
 subroutine math_disturbOri (phi1, PHI, phi2, scatter)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) phi1, PHI, phi2, scatter
 real(pReal) orot(3,3), srot(3,3), p1(2), P(2), p2(2), rot(3,3)
 real(pReal) gscatter,scale,x,y,s,z,arg,angle,rand,gauss

 p1(1)=0
 P(1)=0
 p2(1)=0

! Helming uses different distribution with Bessel functions
! therefore the gauss scatter width has to be scaled differently
 gscatter=0.95*scatter
 scale=cos(gscatter*inRad)

 100 call random_number(x)
 call random_number(y)
 call random_number(s)
 call random_number(z)
 x=x-0.5
 s=s-0.5
 z=z-0.5
 p1(2)=x*gscatter*2.0_pReal
 p2(2)=z*gscatter*2.0_pReal
 arg=scale+y*(1.0-scale)
 P(2)=sign(1.0_pReal,s)*acos(arg)*inDeg
 angle = math_disorient(p1,P,p2)
 call random_number(rand)
 gauss=exp(-1.0*(angle/gscatter)**2)
 if(gauss.LT.rand) then
     goto 100
 end if
! calculate rotation matrix for rotation angles
 srot = math_EulertoR(p1(2),p(2),p2(2))
! calculate rotation matrix for original euler angles
 orot = math_EulertoR(phi1,PHI,phi2)      
! rotate originial orientation matrix
 rot=matmul(srot,orot)
! calculate Euler angles for new rotation matrix
 call math_RtoEuler(rot, phi1,PHI,phi2)
 return
 
 END SUBROUTINE


!********************************************************************
! This routine computes one orientation of a fiber component
!********************************************************************
 subroutine math_fiber(alpha1, alpha2,beta1,beta2,scatter,phi1,PHI,phi2)

 use prec, only: pReal, pInt
 implicit none

 real(pReal) alpha1, alpha2,beta1,beta2,scatter, phi1, PHI, phi2
 real(pReal) orot(3,3), srot(3,3), ac(3), as(3),ori(3,3), rrot(3,3)
 real(pReal) a1r,a2r,b1r,b2r,angle,axis_u,axis_v,axis_w,rand,x,y,z,gscatter,scale,gauss
 integer(pInt) i

! convert angles to radians
 a1r=alpha1*inRad
 a2r=alpha2*inRad
 b1r=beta1*inRad
 b2r=beta2*inRad
! calculate fiber axis in crystal coordinate system
 ac(1)=sin(a1r)*cos(a2r)
 ac(2)=sin(a1r)*sin(a2r)
 ac(3)=cos(a1r)
! calculate fiber axis in sample coordinate system
 as(1)=sin(b1r)*cos(b2r)
 as(2)=sin(b1r)*sin(b2r)
 as(3)=cos(b1r)
! calculate rotation angle between sample and crystal system
 angle=-acos(dot_product(ac, as))
 if(angle.NE.0.0) then
! calculate rotation axis between sample and crystal system
     axis_u=ac(2)*as(3)-ac(3)*as(2)
     axis_v=ac(3)*as(1)-ac(1)*as(3)
     axis_w=ac(1)*as(2)-ac(2)*as(1)
! calculate rotation matrix
     orot = math_RodrigtoR(angle, axis_u, axis_v, axis_w)
 else
     orot = math_I3
 end if

! calculate random rotation angle about fiber axis
 call random_number(rand)
 angle=rand*2.0_pReal*pi
 rrot = math_RodrigtoR(angle, as(1), as(2), as(3))
! find random axis pependicular to fiber axis
 call random_number(x)
 call random_number(y)
 if (as(3).NE.0) then
     z=-(x*as(1)+y*as(2))/as(3)
 else if(as(2).NE.0) then
     z=y
     y=-(x*as(1)+z*as(3))/as(2)
 else if(as(1).NE.0) then
     z=x
     x=-(y*as(2)+z*as(3))/as(1)
 end if
! Helming uses different distribution with Bessel functions
! therefore the gauss scatter width has to be scalled differently
 gscatter=0.95*scatter
 scale=cos(2*gscatter*inRad)
! calculate rotation angle
100 call random_number(rand)
 angle=sign(1.0_pReal,rand)*acos(abs(rand)*scale)*inDeg
 call random_number(rand)
 gauss=exp(-1.0*(angle/gscatter)**2)
 if(gauss.LT.rand) then
     goto 100
 end if
! convert angle to radians
 angle=angle*inRad
 srot = math_RodrigtoR(angle, x, y, z)
 ori=matmul(srot, matmul(rrot, orot))
! calculate Euler angles for new rotation matrix
 call math_RtoEuler(ori, phi1,PHI,phi2)
 
 return
 END SUBROUTINE


 END MODULE math

