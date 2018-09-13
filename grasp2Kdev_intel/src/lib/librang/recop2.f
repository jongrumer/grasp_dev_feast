********************************************************************
*                                                                  *
      SUBROUTINE RECOP2(NS,JA1,JA2,K1,K2,KA,IRE,IAT,REC)
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
      IJ=JLIST(JA1)
      S=DBLE(JJQ1(3,IJ))
      IJ=JLIST(JA2)
      SS=DBLE(JJQ1(3,IJ))
      SS=S*SS
      REC=ONE/SQRT(SS)
      IF(NPEEL.EQ.1 .AND. NS.EQ.-1)RETURN
      IF(NS .EQ. -1) THEN
         NPEELGG = NPEEL
      ELSE
         NPEELGG = NS
      END IF
      IAT=0
      ISKR=NPEELGG-JA2
      IF(ISKR.GT.1) THEN
        CALL DIAGA3(JA2,NPEELGG,2*KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
        IAT=0
      END IF
      IF(JA2.NE.NPEELGG) THEN
        CALL DIAGA5(NPEELGG,JA2,2*KA,IRE,IAT,RE)
        IF(IAT.EQ.0)RETURN
        REC=RE*REC
        IAT=0
      ENDIF
      CALL DIAGA4(JA1,JA2,K1,K2,2*KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      IF(JA1.EQ.1.AND.JA2.EQ.2)RETURN
      IAT=0
      CALL DIAGA1(JA1,K1,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      ISKR=JA2-JA1
      IF(JA1.EQ.1)ISKR=JA2-1-JA1
      IF(ISKR.LE.1)RETURN
      IAT=0
      CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
      REC=RE*REC
      RETURN
      END
