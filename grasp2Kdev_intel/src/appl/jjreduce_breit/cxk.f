************************************************************************
*                                                                      *
      SUBROUTINE CXK (S,IS,KAPS,NU,K,IBR,IEX)
*                                                                      *
*   Computes  the  coefficients of radial integrals in the expansion   *
*   of the effective interaction strength: X(K,IA1,IB1,IA2,IB2).       *
*                                                                      *
*   Input variables:                                                   *
*                                                                      *
*      IS  : Orbital labels                                            *
*      KAPS: Values of 2*kappa                                         *
*      NU  : Order of radial integral                                  *
*      K   : Index of tensor operator                                  *
*      IEX : 1 for direct, 2 for exchange terms                        *
*      IBR : Classifies type of radial integral.There are 4 distinct   *
*            cases:                                                    *
*            IBR = 1 A. All states distinct                            *
*                    B. ((IA .EQ. IB) .AND. (IC .NE. ID)), or          *
*                       ((IA .NE. IB) .AND. (IC .EQ. ID))              *
*                    These give 12 distinct  radial  integrals, with   *
*                    values of K and NU limited only by angular mom-   *
*                    entum and parity                                  *
*            IBR = 2 ((IA .EQ. IC) .AND. (IB .NE. ID)) or              *
*                    ((IA .NE. IC) .AND. (IB .EQ. ID))                 *
*                    This case gives one non-zero integral when K =    *
*                    NU is ODD                                         *
*            IBR = 3 ((IA .EQ. IC) .AND. (IB .EQ. ID)) AND             *
*                    (IA .NE. IB)                                      *
*                    Integrals of magnetic F-type when K = NU is odd   *
*            IBR = 4 ((IA .EQ. ID) .AND. (IB .EQ. IC)) gives 3  mag-   *
*                    netic G-type integrals and  four  H-TYPE  inte-   *
*                    grals                                             *
*                                                                      *
*   Output:                                                            *
*                                                                      *
*      S   : Coefficients S(MU) MU = 1,12                              *
*                                                                      *
*                                                                      *
*   Call(s) to: [LIB92] CRE.                                           *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      DIMENSION IS(4),KAPS(4),S(12)
*
*   1.0  Initialization
*
      DO 1 MU = 1,12
         S(MU) = 0.0D 00
    1 CONTINUE
*
      IA = IS(1)
      IB = IS(2)
      IC = IS(3)
      ID = IS(4)
      KA = KAPS(1)/2
      KB = KAPS(2)/2
      KC = KAPS(3)/2
      KD = KAPS(4)/2
      IF (IEX .NE. 2) GOTO 2
      KK = KD
      IK = ID
      KD = KC
      ID = IC
      KC = KK
      IC = IK
    2 GOTO (3,8,11,12),IBR
      GOTO 17
*
*   2.0  IBR = 1 --- The general case
*
    3 CONTINUE
      IF (NU-K) 7,4,6
*
*   2.1  NU = K .GT. 0
*
    4 CONTINUE
      S(1) = -(KA+KC)*(KD+KB)
      IF (K .EQ. 0) GOTO 16
      D = K*(K+1)
      H = CRE (KA,K,KC)*CRE(KB,K,KD)
      IF (MOD (K,2) .NE. 0) H = -H
      S(1) = S(1)*H/D
      DO 5 MU = 2,4
         S(MU) = S(1)
    5 CONTINUE
      RETURN
*
*   2.2  NU = K+1
*
    6 CONTINUE
      DK1 = KC-KA
      DK2 = KD-KB
      FK = K
      GK = K+1
      G1 = DK1-GK
      G2 = DK1+GK
      G3 = DK2-GK
      G4 = DK2+GK
      KK = K+K+1
      H = CRE (KA,K,KC)*CRE(KB,K,KD)
      IF (MOD (K,2) .NE. 0) H = -H
      A = H*FK/GK/DBLE (KK*(KK+2))
      S(1) = A*G1*G3
      S(2) = A*G2*G4
      S(3) = A*G1*G4
      S(4) = A*G2*G3
      RETURN
*
*   2.2  NU = K-1
*
    7 CONTINUE
      DK1 = KC-KA
      DK2 = KD-KB
      FK = K
      GK = K+1
      F1 = DK1-FK
      F2 = DK1+FK
      F3 = DK2-FK
      F4 = DK2+FK
      G1 = DK1-GK
      G2 = DK1+GK
      G3 = DK2-GK
      G4 = DK2+GK
      KK = K+K+1
      H = CRE (KA,K,KC)*CRE(KB,K,KD)
      IF (MOD (K,2) .NE. 0) H = -H
      A = H*GK/FK/DBLE (KK*(KK-2))
      S(1) = A*F2*F4
      S(2) = A*F1*F3
      S(3) = A*F2*F3
      S(4) = A*F1*F4
      B = H/DBLE (KK*KK)
      S(5) = B*F2*G3
      S(6) = B*F4*G1
      S(7) = B*F1*G4
      S(8) = B*F3*G2
      S(9) = B*F2*G4
      S(10) = B*F3*G1
      S(11) = B*F1*G3
      S(12) = B*F4*G2
      RETURN
*
*   3.0  IBR = 2  Degenerate case: only one non-zero R-integral
*
    8 CONTINUE
      IF ((IA .EQ. IC) .AND. (IB .NE. ID)) GOTO 10
      IF ((IA .NE. IC) .AND. (IB .EQ. ID)) GOTO 9
      GOTO 17
*
    9 IK = IB
      IB = IA
      IA = IK
      IK = ID
      ID = IC
      IC = IK
*
      KK = KB
      KB = KA
      KA = KK
      KK = KD
      KD = KC
      KC = KK
*
   10 IF (MOD (K,2) .NE. 1) RETURN
      DK = K*(K+1)
      H = CRE (KA,K,KC)*CRE(KB,K,KD)/DK
      S(1) = H*DBLE (4*KA*(KB+KD))
      RETURN
*
*   4.0  IBR = 3. Direct magnetic F-integrals
*
   11 CONTINUE
      IF ((IA .NE. IC) .OR. (IB .NE. ID)) GOTO 17
      IF (MOD (K,2) .NE. 1) RETURN
      DK = K*(K+1)
      H = CRE(KA,K,KA)*CRE(KB,K,KB)/DK
      S(1) = H*DBLE (16*KA*KB)
      RETURN
*
*   5.0   IBR = 4. Exchange magnetic G- and H-integrals
*
   12 CONTINUE
      IF ((IA .NE. ID) .OR. (IB .NE. IC) )GOTO 17
      IF (NU-K) 15,13,14
*
*   5.1  NU = K
*
   13 CONTINUE
      S(1) = DBLE (KA+KB)*CRE(KA,K,KB)
      IP = ABS (KA)-ABS (KB)+K+1
      S(1) = S(1)*S(1)/DBLE(K*(K+1))
      IF (MOD (IP,2) .NE. 0) S(1) = -S(1)
      S(3) = S(1)
      S(2) = S(1)+S(1)
      RETURN
*
*   5.2  NU = K+1
*
   14 CONTINUE
      DK = KB-KA
      GK = K+1
      FK = K
      G1 = DK+GK
      G2 = DK-GK
      KK = K+K+1
      H = CRE (KA,K,KB)**2
      IF (KA*KB .LT. 0) H = -H
      A = H*FK/GK/DBLE(KK*(KK+2))
      S(1) = -A*G1*G1
      S(2) = -2.0D 00*A*G1*G2
      S(3) = -A*G2*G2
      RETURN
*
*   5.3  NU = K-1
*
   15 CONTINUE
      DK = KB-KA
      FK = K
      GK = K+1
      F1 = DK+FK
      F2 = DK-FK
      G1 = DK+GK
      G2 = DK-GK
      KK = K+K+1
      H = CRE (KA,K,KB)**2
      IF (KA*KB .LT. 0) H = -H
      A = H*GK/FK/DBLE (KK*(KK-2))
      S(1) = -A*F2*F2
      S(2) = -2.0D 00*A*F1*F2
      S(3) = -A*F1*F1
      B = H/DBLE (KK*KK)
      B = B+B
      S(4) = -B*F1*G2
*     S(5) = S(4)
      S(5) = -B*F2*G1
      S(6) = -B*F1*G1
      S(7) = -B*F2*G2
      RETURN
*
*   6.0  Special cases and errors
*
*   Illegal zero value of K in Type 1
*
   16 WRITE (*,300) IS(1),IS(2),IS(3),IS(4),NU,IBR,IEX
      STOP
*
*   Illegal combination of states in Type 3 or 4
*
   17 WRITE (*,301) IBR,IS(1),IS(2),IS(3),IS(4),NU,K,IEX
      STOP
*
  300 FORMAT ('CXK: Illegal value K = 0 -'
     :       /1X,4I3,2X,I3,2X,2I2)
  301 FORMAT ('CXK: Type ',I2,'-'
     :       /1X,I2,3X,4I3,2X,2I3,2X,I2)
*
      END
