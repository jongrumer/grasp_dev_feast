********************************************************************
*                                                                  *
      SUBROUTINE RECOP1(NS,JA1,KA,IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 06  --------------  *
*                                                                  *
*     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3                     *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :       /M0/JJC1(NNNW),JJC2(NNNW)
     :       /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      IAT=1
      IJ1=JLIST(JA1)
      S=DBLE(JJQ1(3,IJ1))
      REC=ONE/SQRT(S)
      IF(NPEEL.EQ.1 .AND. NS.EQ.-1)RETURN
      IF(NS .EQ. -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IAT=0
      IF(IRE.NE.0) THEN
        IF(KA.EQ.0) THEN
          IAT=1
          RETURN
        END IF
      END IF
      IAT=1
      IF(NPEELGG.EQ.1) RETURN
      IAT=0
      IF(NPEELGG.NE.2) THEN
        CALL DIAGA5(NPEELGG,JA1,2*KA,IRE,IAT,RE)
        REC=RE*REC
        IF(IAT.EQ.0) RETURN
        IF(JA1.EQ.NPEELGG) RETURN
        IAT=0
      END IF
      CALL DIAGA1(JA1,2*KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      IF(NPEELGG.EQ.2) RETURN
      ISKR=NPEELGG-JA1
      IF(JA1.EQ.1)ISKR=NPEELGG-1-JA1
      IF(ISKR.LE.1)RETURN
      IAT=0
      CALL DIAGA3(JA1,NPEELGG,2*KA,IRE,IAT,RE)
      REC=RE*REC
      RETURN
      END
