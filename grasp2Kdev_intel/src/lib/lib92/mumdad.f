************************************************************************
*                                                                      *
      SUBROUTINE MUMDAD (IS,KAPS,X)
*                                                                      *
*   Evaluate the product of 4 CFPs.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: CFP.                                          *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
*
      PARAMETER (EPS = 1.0D-10)
*
      DIMENSION IS(2,2),KAPS(2,2)
*
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /L1/JBQ1(3,NNNW),JBQ2(3,NNNW),JTQ1(3),JTQ2(3)
*
      X = 1.0D 00
*
*   First index
*
      LOCK = KAPS(1,1)
      IF (ABS (LOCK) .EQ. 2) GOTO 4
      II = IS(1,1)
      NEL = NQ1(II)
      IVP = JBQ1(1,II)
      IWP = JBQ1(2,II)
      IJP = JBQ1(3,II)-1
*
*   IA1 .NE. IB1 and IA2 .NE. IB2; use JJQ array.
*
      IF (IS(1,1) .EQ. IS(2,1)) GOTO 1
      IVD = JJQ1(1,II)
      IWD = JJQ1(2,II)
      IJD = JJQ1(3,II)-1
      GOTO 2
*
*   IA1 .EQ. IB1 or IA2 .EQ. IB2; JTQ array needed.
*
    1 NEL = NEL-1
      IVD = JTQ1(1)
      IWD = JTQ1(2)
      IJD = JTQ1(3)-1
    2 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
*
    4 LOCK = KAPS(2,1)
      IF (IABS (LOCK) .EQ. 2) GOTO 8
      II = IS(2,1)
      NEL = NQ1(II)
      IVD = JJQ1(1,II)
      IWD = JJQ1(2,II)
      IJD = JJQ1(3,II)-1
      IF (IS(1,1) .EQ. IS(2,1)) GOTO 5
      IVP = JBQ1(1,II)
      IWP = JBQ1(2,II)
      IJP = JBQ1(3,II)-1
      GOTO 6
    5 IVP = JTQ1(1)
      IWP = JTQ1(2)
      IJP = JTQ1(3)-1
    6 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
    8 CONTINUE
*
*   Second index
*
      LOCK = KAPS(1,2)
      IF (ABS (LOCK) .EQ. 2) GOTO 12
      II = IS(1,2)
      NEL = NQ2(II)
      IVP = JBQ2(1,II)
      IWP = JBQ2(2,II)
      IJP = JBQ2(3,II)-1
*
*   IA1 .NE. IB1 and IA2 .NE. IB2; use JJQ array.
*
      IF (IS(1,2) .EQ. IS(2,2)) GOTO 9
      IVD = JJQ2(1,II)
      IWD = JJQ2(2,II)
      IJD = JJQ2(3,II)-1
      GOTO 10
*
*   IA1 .EQ. IB1 or IA2 .EQ. IB2; JTQ array needed.
*
    9 NEL = NEL-1
      IVD = JTQ2(1)
      IWD = JTQ2(2)
      IJD = JTQ2(3)-1
   10 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
*
   12 LOCK = KAPS(2,2)
      IF (ABS (LOCK) .EQ. 2) GOTO 16
      II = IS(2,2)
      NEL = NQ2(II)
      IVD = JJQ2(1,II)
      IWD = JJQ2(2,II)
      IJD = JJQ2(3,II)-1
      IF (IS(1,2) .EQ. IS(2,2)) GOTO 13
      IVP = JBQ2(1,II)
      IWP = JBQ2(2,II)
      IJP = JBQ2(3,II)-1
      GOTO 14
   13 IVP = JTQ2(1)
      IWP = JTQ2(2)
      IJP = JTQ2(3)-1
   14 CALL CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C)
      IF (IBUG2 .NE. 0)
     :   WRITE (99,300) LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,C
      IF (ABS (C) .LT. EPS) GOTO 17
      X = X*C
   16 CONTINUE
      RETURN
*
   17 X = 0.0D 00
      RETURN
*
  300 FORMAT ('MUMDAD: CFP ',I3,I4,I7,2I4,I7,2I4,1P,1D19.12)
*
      END
