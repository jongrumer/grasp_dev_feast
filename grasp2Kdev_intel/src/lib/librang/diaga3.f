********************************************************************
*                                                                  *
      SUBROUTINE DIAGA3(JA1,JA2,KA,IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 03  --------------  *
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
      AA=ONE
      I=JA1+1
      IF(JA1.EQ.1)I=I+1
      IF(I.GE.JA2)GO TO 4
    1 LL1=JLIST(I)
      JI=JJQ1(3,LL1)-1
      KK2=I-2
      ITI=JJC1(KK2)-1
      ITIS=JJC2(KK2)-1
      KK1=I-1
      ITI1=JJC1(KK1)-1
      ITI1S=JJC2(KK1)-1
      IF(IRE.NE.0)GO TO 2
      IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S).EQ.0)RETURN
      GO TO 3
    2 CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
      A3=A3*SQRT(DBLE((ITI+1)*(ITI1S+1)))
      IFAZ=KA+JI+ITI+ITI1S
      IF((IFAZ/4)*4.NE.IFAZ)A3=-A3
      AA=AA*A3
    3 I=I+1
      IF(I.EQ.JA2)GO TO 4
      GO TO 1
    4 REC=AA
      IAT=1
      RETURN
      END
