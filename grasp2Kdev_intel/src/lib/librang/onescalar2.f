********************************************************************
*                                                                  *
      SUBROUTINE ONESCALAR2(JJA,JJB,JA,JB,COEFF)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
*                                                                  *
*     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 - 1   *
*                                                   N'2 = N2 + 1   *
*                                                                  *
*     SUBROUTINE CALLED: COULOM,GG1222,ITREXG,IXJTIK,PERKO2,       *
*                        RECO,RECO2,SIXJ,SPEAK                     *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/M1/NQ1(NNNW),NQ2(NNNW)
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
      COEFF=ZERO
      IF(JA.EQ.JB) RETURN
      IF(JA.LT.JB) THEN
        JAA=JA
        JBB=JB
      ELSE
        JAA=JB
        JBB=JA
      END IF
      CALL RECOONESCALAR(-1,JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT.EQ.0)RETURN
      QM1=HALF
      QM2=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      IF(ID1(3).NE.ID2(3)) RETURN
      CALL RECO2(JAA,JBB,ID2(3),0,IAT,REC)
      IF(IAT.EQ.0)RETURN
      CALL GG12(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,QM1,QM2,WW)
      IF(DABS(WW).GT.EPS) THEN
         CALL RECO2(JAA,JBB,ID2(3),1,IAT,REC)
         COEFF=WW*REC*DSQRT(DBLE(ID1(3)+1))
         NN=0
         IB1=JBB-1
         DO 1 II=JAA,IB1
           IN=JLIST(II)
           NN=NQ1(IN)+NN
    1    CONTINUE
         IF((NN/2)*2.EQ.NN) COEFF=-COEFF
CGG         IF(JA.GT.JB) COEFF=-COEFF
         COEFF=-COEFF
C
C     TRANSFORM FANO & RACAH PHASE CONVENTION
C     TO CONDON & SHORTLEY PHASE CONVENTION
C
        IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
        IF((IFAZ/4)*4.NE.IFAZ)COEFF=-COEFF
      END IF
      RETURN
      END
