********************************************************************
*                                                                  *
      SUBROUTINE RECO4(JA1,JA2,JA3,JA4,K1,K2,K3,K4,KA,
     *IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 09  --------------  *
*                                                                  *
*     SUBROUTINE CALLED:  DIAGA1,DIAGA2,DIAGA3,DIAGA4              *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :       /M0/JJC1(NNNW),JJC2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IJ3=JLIST(JA3)
      IJ4=JLIST(JA4)
      S1=JJQ1(3,IJ1)
      S2=JJQ1(3,IJ2)
      S3=JJQ1(3,IJ3)
      S4=JJQ1(3,IJ4)
      S=S1*S2*S3*S4
      REC=1.0/SQRT(S)
      IA4=JJQ1(3,IJ4)-1
      IB4=JJQ2(3,IJ4)-1
      REC=REC*SQRT(DBLE(IA4+1))/
     *SQRT(DBLE((K4+1)*(IB4+1)))
*
      ISKR=JA3-JA2
      IF(ISKR.LE.1)GO TO 1
      IAT=0
      CALL DIAGA3(JA2,JA3,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
    1 ISKR=JA4-JA3
      IF(ISKR.LE.1)GO TO 2
      IAT=0
      CALL DIAGA3(JA3,JA4,K4,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
    2 IAT=0
      CALL DIAGA2(JA1,JA4,K4,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      IAT=0
      CALL DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      IAT=0
      CALL DIAGA4(JA2,JA3,KA,K3,K4,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      IF(JA1.EQ.1.AND.JA2.EQ.2)RETURN
*
      IAT=0
      CALL DIAGA1(JA1,K1,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
*
      ISKR=JA2-JA1
      IF(JA1.EQ.1)ISKR=JA2-1-JA1
      IF(ISKR.LE.1)RETURN
      IAT=0
      CALL DIAGA3(JA1,JA2,K1,IRE,IAT,RE)
      REC=RE*REC
      RETURN
      END
