********************************************************************
*                                                                  *
      SUBROUTINE RECO2(JA1,JA2,KA,IRE,IAT,REC)
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
      IAT=0
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      S=DBLE(JJQ1(3,IJ1))
      SS=DBLE(JJQ1(3,IJ2))
      SS=S*SS
      REC=ONE/SQRT(SS)
      IF(IRE.EQ.0)GO TO 1
      IF(KA.NE.0)GO TO 1
      IAT=1
      RETURN
*
    1 IA1=JJQ1(3,IJ1)-1
      IA2=JJQ1(3,IJ2)-1
      IB1=JJQ2(3,IJ1)-1
      IB2=JJQ2(3,IJ2)-1
      IAT=0
      CALL DIAGA2(JA1,JA2,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC*SQRT(DBLE(IA2+1))/
     *SQRT(DBLE((KA+1)*(IB2+1)))
      IF(JA1.EQ.1.AND.JA2.EQ.2)RETURN
*
      IAT=0
      CALL DIAGA1(JA1,KA,IRE,IAT,RE)
      IF(IAT.EQ.0)RETURN
      REC=RE*REC
      ISKR=JA2-JA1
      IF(JA1.EQ.1)ISKR=JA2-1-JA1
      IF(ISKR.LE.1)RETURN
*
      IAT=0
      CALL DIAGA3(JA1,JA2,KA,IRE,IAT,RE)
      REC=RE*REC
      RETURN
      END
