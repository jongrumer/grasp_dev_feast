********************************************************************
*                                                                  *
      SUBROUTINE DIAGA5(NPEELGG,JA1,KA,IRE,IAT,REC)
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
      ITI1=JJC1(NPEELGG-1)-1
      ITI1S=JJC2(NPEELGG-1)-1
      IJ1=JLIST(NPEELGG)
      IF(JA1.EQ.NPEELGG) THEN
        ITI=JJQ1(3,IJ1)-1
        ITIS=JJQ2(3,IJ1)-1
        JI=JJC1(NPEELGG-2)-1
      ELSE
        JI=JJQ1(3,IJ1)-1
        ITI=JJC1(NPEELGG-2)-1
        ITIS=JJC2(NPEELGG-2)-1
      END IF
      IF(IRE.EQ.0) THEN
        IF(IXJTIK(KA,ITIS,ITI,JI,ITI1,ITI1S).NE.0) IAT=1
      ELSE
        CALL SIXJ(KA,ITIS,ITI,JI,ITI1,ITI1S,0,A3)
        REC=A3*SQRT(DBLE((ITI+1)*(ITI1S+1)))
        IF(MOD(KA+JI+ITIS+ITI1,4).NE.0)REC=-REC
        IAT=1
        IF(JA1.EQ.NPEELGG)RETURN
        IF(MOD(ITI+ITIS-ITI1S-ITI1+2*JI,4).NE.0)REC=-REC
      END IF
      RETURN
      END
