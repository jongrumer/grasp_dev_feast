************************************************************************
*                                                                      *
      SUBROUTINE SETQNA (JA,JB)
*                                                                      *
*   This generates the  arrays  defining  the quantum numbers of the   *
*   states involved in the  matrix  element  linking  configurations   *
*   labelled by JA, JB.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: ICHOP, IQ, JCUP, JQS.                         *
*                                                                      *
*                                           Last update: 30 Oct 1987   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/M0/JJC1(NNNW),JJC2(NNNW)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
*
*   List parameters defining all shells in both configurations, whether
*   participating or not
*
      DO 2 J = 1,NW
         NQ1(J) = IQ (J,JA)
         NQ2(J) = IQ (J,JB)
         DO 1 K = 1,3
            JJQ1(K,J) = JQS (K,J,JA)
            JJQ2(K,J) = JQS (K,J,JB)
    1    CONTINUE
    2 CONTINUE
*
*   Define coupling schemes: set JLIST array to define those shells
*   which are open in either configuration, and KLIST array to locate
*   the rest. Exclude shells which are empty in both configurations
*
      NPEEL = 0
      NCORE = 0
      DO 3 J = 1,NW
         IF ((ICHOP (J,JA) .NE. -1) .OR. (ICHOP (J,JB) .NE. -1)) THEN
            IF ((ICHOP (J,JA) .EQ. 1) .AND. (ICHOP (J,JB) .EQ. 1))
     :         THEN
               NCORE = NCORE+1
               KLIST(NCORE) = J
            ELSE
               NPEEL = NPEEL+1
               JLIST(NPEEL) = J
            ENDIF
         ENDIF
    3 CONTINUE
*
*   Return if not more than one shell is open
*
      IF (NPEEL .LE. 1) RETURN
*
*   Set arrays of coupling angular momenta interpolating closed
*   shells where necessary. Left hand side first ...
*
      JCNT = 1
      JCNTOP = 0
      JW1 = JLIST(1)
      JW2 = JLIST(2)
      IF (ICHOP (JW1,JA) .NE. 0) THEN
         JJC1(1) = JQS (3,JW2,JA)
         IF (ICHOP (JW2,JA) .EQ. 0) JCNTOP = 1
      ELSE
         JCNTOP = 1
         IF (ICHOP (JW2,JA) .EQ. 0) THEN
            JJC1(1) = JCUP (JCNT,JA)
            JCNT = JCNT+1
         ELSE
            JJC1(1) = JQS (3,JW1,JA)
         ENDIF
      ENDIF
*
      DO 4 J = 3,NPEEL
         JW = JLIST(J)
         IF (ICHOP (JW,JA) .NE. 0) THEN
            JJC1(J-1) = JJC1(J-2)
         ELSE
            IF (JCNTOP .NE. 0) THEN
               JJC1(J-1) = JCUP (JCNT,JA)
               JCNT = JCNT+1
            ELSE
               JJC1(J-1) = JQS (3,JW,JA)
            ENDIF
            JCNTOP = JCNTOP+1
         ENDIF
    4 CONTINUE
*
*   ... and repeat for right hand side
*
      JCNT = 1
      JCNTOP = 0
      JW1 = JLIST(1)
      JW2 = JLIST(2)
      IF (ICHOP (JW1,JB) .NE. 0) THEN
         JJC2(1) = JQS (3,JW2,JB)
         IF (ICHOP (JW2,JB) .EQ. 0) JCNTOP = 1
      ELSE
         JCNTOP = 1
         IF (ICHOP (JW2,JB) .EQ. 0) THEN
            JJC2(1) = JCUP (JCNT,JB)
            JCNT = JCNT+1
         ELSE
            JJC2(1) = JQS (3,JW1,JB)
         ENDIF
      ENDIF
*
      DO 5 J = 3,NPEEL
         JW = JLIST(J)
         IF (ICHOP (JW,JB) .NE. 0) THEN
            JJC2(J-1) = JJC2(J-2)
         ELSE
            IF (JCNTOP .NE. 0) THEN
               JJC2(J-1) = JCUP (JCNT,JB)
               JCNT = JCNT+1
            ELSE
               JJC2(J-1) = JQS (3,JW,JB)
            ENDIF
            JCNTOP = JCNTOP+1
         ENDIF
    5 CONTINUE
*
      RETURN
      END
