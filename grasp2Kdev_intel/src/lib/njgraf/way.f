************************************************************************
*                                                                      *
      SUBROUTINE WAY (L,KA,KB,ICH,NB)
*                                                                      *
*   Tests  one  step  forward  if  the way is free. First and second   *
*   arguments are interchanged or not according to ICH = -1, or +1.    *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INTEGER ARROW
      LOGICAL TABS
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
*
      K1 = J23(L,KA)
      K2 = J23(L,KB)
      NB = IAL(K1)+IAL(K2)-1
      IF (NB) 3,2,8
    2 NB1 = IAL(K1)-IAL(K2)
      IF (NB1) 9,8,8
    3 CALL OTHERJ (L,K1,L1,LC1,LA)
      CALL OTHERJ (L,K2,L2,LC2,LB)
      CALL NEIBOR (LC1,I1,I2)
      CALL NEIBOR (LC2,I3,I4)
      JI1 = J23(L1,I1)
      JI2 = J23(L1,I2)
      JI3 = J23(L2,I3)
      JI4 = J23(L2,I4)
      IA = IAL(JI1)+IAL(JI2)
      IB = IAL(JI3)+IAL(JI4)
      NBP = IB+IA+1
      NBM = IB-IA
      GOTO (8,4,5,4,6),NBP
    4 IF (NBM) 9,8,8
    5 IF (NBM) 9,6,8
    6 IF ((JI3 .EQ. IF1) .OR. (JI3 .EQ. IF2) .OR.
     :    (JI4 .EQ. IF1) .OR. (JI4 .EQ. IF2)) GOTO 9
    8 ICH = 1
      GOTO 10
    9 ICH = -1
   10 RETURN
*
      END
