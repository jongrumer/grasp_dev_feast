********************************************************************
*                                                                  *
      SUBROUTINE DIAGA1(JA1,KA,IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 01  --------------  *
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
      IA1=JJQ1(3,IJ1)-1
      IB1=JJQ2(3,IJ1)-1
      IF(JA1.LE.2)GO TO 1
      K1=JA1-2
      J1=JJC1(K1)-1
      K1=JA1-1
      IT1=JJC1(K1)-1
      IT1S=JJC2(K1)-1
      GO TO 2
    1 LL1=JLIST(1)
      IF(JA1.EQ.1)LL1=JLIST(2)
      J1=JJQ1(3,LL1)-1
      IT1=JJC1(1)-1
      IT1S=JJC2(1)-1
    2 IF(IRE.NE.0)GO TO 3
      IF(IXJTIK(KA,IB1,IA1,J1,IT1,IT1S).EQ.0)RETURN
      IAT=1
      RETURN
    3 CALL SIXJ(KA,IB1,IA1,J1,IT1,IT1S,0,A1)
      A1=A1*SQRT(DBLE((IA1+1)*(IT1S+1)))
      IFAZ=J1+IT1+IB1+KA
      IF((IFAZ/4)*4.NE.IFAZ)A1=-A1
      REC=A1
      IAT=1
      IF(JA1.NE.1)RETURN
      IFAZ=IA1+IB1+2*J1-IT1-IT1S
      IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
      RETURN
      END
