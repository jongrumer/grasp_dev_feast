********************************************************************
*                                                                  *
      SUBROUTINE RECOONESCALAR(NS,JA1,JA2,JA3,JA4,KA,IAT)
*                                                                  *
*     -------------  SECTION REC    SUBPROGRAM 05  --------------  *
*                                                                  *
*     NO SUBROUTINE CALLED                                         *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
      COMMON/M1/NQ1(NNNW),NQ2(NNNW)
     :       /M0/JJC1(NNNW),JJC2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
      IAT=1
      IF(NPEEL.EQ.1 .AND. NS.EQ.-1)RETURN
      IF(NS .EQ. -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IF(JA1.EQ.1.AND.JA2.EQ.2) GO TO 1
      IF(KA.NE.0)GO TO 5
*
*  CASES WHEN :          KA = 0
*                  OR    JA1 = JA2
*                  OR    JA1 = 1    JA2 = 2
*
    1 DO 3 I=1,NPEELGG
        IJ=JLIST(I)
        IF(I.GE.NPEELGG-1)GO TO 4
        IF(JJC1(I).NE.JJC2(I))IAT=0
    4   IF(KA.EQ.0)GO TO 9
        IF(I.EQ.JA1)GO TO 3
        IF(I.EQ.JA2)GO TO 3
    9   CONTINUE
        DO 2 J=1,3
          IF(JJQ1(J,IJ).NE.JJQ2(J,IJ))IAT=0
    2   CONTINUE
    3 CONTINUE
      RETURN
*
*  OTHER CASES
*
    5 CONTINUE
      DO 6 I=1,NPEELGG
        IJ=JLIST(I)
        IF(I.GE.NPEELGG-1)GO TO 7
        IA1=JA1-1
        IA2=JA2-1
        IF(JA1.EQ.1)IA1=JA1
        IF(I.GE.IA1.AND.I.LT.IA2)GO TO 7
        IF(JJC1(I).NE.JJC2(I))IAT=0
    7   IF(I.EQ.JA1)GO TO 6
        IF(I.EQ.JA2)GO TO 6
        IF((KA.EQ.2).AND.(I.EQ.JA3))GO TO 6
        IF((KA.EQ.3).AND.(I.EQ.JA3))GO TO 6
        IF((KA.EQ.3).AND.(I.EQ.JA4))GO TO 6
        DO 8 J=1,3
          IF(JJQ1(J,IJ).NE.JJQ2(J,IJ))IAT=0
    8   CONTINUE
    6 CONTINUE
      RETURN
      END
