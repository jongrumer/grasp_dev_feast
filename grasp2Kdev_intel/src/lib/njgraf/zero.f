************************************************************************
*                                                                      *
      SUBROUTINE ZERO (J,JZ,FAIL)
*                                                                      *
*   Suppresses  one  line  and  two  nodes of the unstructured graph   *
*   introduces  zeros in the triads  J23. As a consequence the other   *
*   two arguments of the triad are put equal. If there was already a   *
*   zero in the triad which is changed, it is a special case.          *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW
      LOGICAL FAIL,TABS,FREE,SUMVAR,CUT,NOCUT
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10,MZERO = 20)
*
      COMMON/ZER/NZERO,JZERO(MZERO)
     :      /CUTDIG/CUT
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'ZERO  '/
*
      NOCUT = .FALSE.
      NZERO = 0
*
      IF (J .GE. 1) THEN
         CALL OTHERJ (0,JZ,LIN,LC,K1)
         I = NZERO
         GOTO 8
      ENDIF
*
      DO 11 I = 1,M
         IF ((J1(I) .NE. 1) .OR. FREE(I) .OR. (IAL(I).LE.1)) GOTO 11
         NZERO = NZERO+1
         IF (NZERO .GT. MZERO) THEN
            WRITE (*,300) NZERO,MZERO
            STOP
         ENDIF
         JZERO(NZERO) = I
   11 CONTINUE
*
      NOCUT = .TRUE.
      M = M+1
      J1(M) = 1
      SUMVAR(M) = .FALSE.
      FREE(M) = .FALSE.
      IF (NZERO .EQ. 0) GOTO 7
      CALL PRINTJ (NAME,1)
      I = 0
    1 I = I+1
      JZ = JZERO(I)
      J = 0
   13 J = J+1
      LIN = LINE(JZ,J)
      IF (TABS(LIN)) GOTO 2
      LC = LCOL(JZ,J)
    8 CALL NEIBOR (LC,L1,L2)
      JJ1 = J23(LIN,L1)
      JJ2 = J23(LIN,L2)
*
      IF (JJ1 .EQ. JJ2) THEN
         J6C = J6C+1
         J6(J6C) = JJ1
         LO1=LIN
         LO2=LIN
         LCO1=L1
         LCO2=L2
         GOTO 10
      ENDIF
*
      CALL DELTA (JJ1,JJ2,FAIL)
      IF (FAIL) GOTO 7
*
      IF ((J1(JJ1) .NE. 1) .AND. (J1(JJ2) .NE. 1)) GOTO 15
      IF (J1(JJ1) .LT. J1(JJ2)) GOTO 15
      IF (J1(JJ1) .GT. J1(JJ2)) GOTO 19
*
      IF (NZERO .NE. 0) THEN
         DO 17 JJX = I,NZERO
            JJZ = JZERO(JJX)
            IF (JJ1 .EQ. JJZ) GOTO 15
            IF (JJ2 .EQ. JJZ) GOTO 19
   17    CONTINUE
      ENDIF
*
      GOTO 15
*
   19 JJZ = JJ2
      JJ2 = JJ1
      JJ1 = JJZ
*
   15 CALL OTHERJ (LIN,JJ1,LO1,LCO1,K1)
      CALL OTHERJ (LIN,JJ2,LO2,LCO2,K2)
      J9C = J9C+1
      J9(J9C) = JJ1
      J23(LO2,LCO2) = JJ1
      LINE(JJ1,K1) = LO2
      LCOL(JJ1,K1) = LCO2
*
   10 IF     (ARROW(LIN,L1) .LT. ARROW(LIN,L2)) THEN
         CALL PHASE2 (JJ1)
      ELSEIF (ARROW(LIN,L1) .EQ. ARROW(LIN,L2)) THEN
         ARROW(LO1,LCO1) = 1
         ARROW(LO2,LCO2) = -1
      ENDIF
*
      TABS(LIN) = .TRUE.
      NBTR = NBTR-1
      IF (NBTR .EQ. 0) GOTO 7
      IF (LO1 .EQ. LO2) THEN
         L = 6-LCO1-LCO2
         JT = J23(LO1,L)
         IF ((J1(JT) .EQ. 1) .AND. (.NOT.FREE(JT))) GOTO 2
         CALL DELTA (JT,M,FAIL)
         IF (FAIL) GOTO 7
         NZERO = NZERO+1
         JZERO(NZERO) = JT
      ENDIF
    2 IF (J .EQ. 1) GOTO 13
*
      IF (NBTR .NE. 0) THEN
         IF (I .LT. NZERO) GOTO 1
      ENDIF
*
    7 CALL PRINTJ (NAME,4)
      IF (NOCUT) CUT = .FALSE.
*
      RETURN
*
  300 FORMAT (' Dimension error in ZERO. NZERO = ',I5,' MZERO = ',I5)
*
      END
