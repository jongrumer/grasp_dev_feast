********************************************************************
*                                                                  *
      SUBROUTINE DIAGA4(JA1,JA2,K1,K2,KA,IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 04  --------------  *
*                                                                  *
*     SUBROUTINE CALLED:  NINE                                     *
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
      IA1=JJQ1(3,IJ1)-1
      IA2=JJQ1(3,IJ2)-1
      IB1=JJQ2(3,IJ1)-1
      IB2=JJQ2(3,IJ2)-1
      IF(JA1.EQ.1.AND.JA2.EQ.2)GO TO 1
      N1=JA2-1
      J2=JJC1(N1)-1
      J2S=JJC2(N1)-1
      N2=JA2-2
      IT2=JJC1(N2)-1
      IT2S=JJC2(N2)-1
      GO TO 2
    1 IT2=IA1
      IT2S=IB1
      J2=JJC1(1)-1
      J2S=JJC2(1)-1
    2 IF(IRE.NE.0)GO TO 3
      CALL NINE(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,1,IAT,A2)
      RETURN
    3 CALL NINE(IT2S,K1,IT2,IB2,K2,IA2,J2S,KA,J2,0,IAT,A2)
      REC=A2*SQRT(DBLE((IT2+1)*(KA+1)*(IA2+1)*(J2S+1)))
      RETURN
      END
