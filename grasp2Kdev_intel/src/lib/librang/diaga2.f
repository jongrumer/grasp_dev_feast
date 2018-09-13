********************************************************************
*                                                                  *
      SUBROUTINE DIAGA2(JA1,JA2,KA,IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 02  --------------  *
*                                                                  *
*     SUBROUTINE CALLED:  IXJTIK, SIXJ                             *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :       /M0/JJC1(NNNW),JJC2(NNNW)
     :       /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      REC = ZERO
      IJ1=JLIST(JA1)
      IJ2=JLIST(JA2)
      IA1=JJQ1(3,IJ1)-1
      IA2=JJQ1(3,IJ2)-1
      IB1=JJQ2(3,IJ1)-1
      IB2=JJQ2(3,IJ2)-1
      IF(JA1.EQ.1.AND.JA2.EQ.2)GO TO 1
      N1=JA2-1
      J2=JJC1(N1)-1
      N2=JA2-2
      IT2=JJC1(N2)-1
      IT2S=JJC2(N2)-1
      GO TO 2
    1 IT2=IA1
      IT2S=IB1
      J2=JJC1(1)-1
    2 IF(IRE.NE.0)GO TO 3
      IF(IXJTIK(KA,IB2,IA2,J2,IT2,IT2S).EQ.0)RETURN
      IAT=1
      RETURN
    3 CALL SIXJ(KA,IB2,IA2,J2,IT2,IT2S,0,A2)
      REC=A2*SQRT(DBLE((IB2+1)*(IT2+1)))
      IFAZ=J2+IT2S+IA2+KA
      IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
      IAT=1
      RETURN
      END
