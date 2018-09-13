********************************************************************
*                                                                  *
      SUBROUTINE ONEPARTICLEJJ1(NS,KA,JJA,JJB,JA,JB,COEFF)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF ONE PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                  N'2 = N2        *
*                                                                  *
*      SUBROUTINE CALLED:                                          *
*                                                                  *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
C
C     THE CASE 11   + -
C
      COEFF=ZERO
      IF(JA.EQ.JB) THEN
        IF(JJA.NE.JJB) THEN
          CALL RECOP00(NS,JA,JA,KA,IAT)
          IF(IAT.EQ.0)RETURN
        END IF
        CALL RECOP1(NS,JA,KA,0,IAT,REC)
        IF(IAT.EQ.0)RETURN
        CALL PERKO2(JA,JA,JA,JA,1)
        QM1=HALF
        QM2=-HALF
        CALL WJ1(IK1,BK1,ID1,BD1,KA,QM1,QM2,WJ)
        IF(DABS(WJ).GT.EPS) THEN
           CALL RECOP1(NS,JA,KA,1,IAT,REC)
           COEFF = WJ*REC*DSQRT(DBLE(ID1(3)+1))
     :           / DSQRT(DBLE(2*KA+1))
           COEFF=-COEFF
        END IF
      END IF
      RETURN
      END
