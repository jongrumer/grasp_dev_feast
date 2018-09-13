********************************************************************
*                                                                  *
      SUBROUTINE WAP1(IK,BK,ID,BD,K1,BK2,QM1,QM2,QM3,WA)
*                                                                  *
*   ---------------  SECTION SQ    SUBPROGRAM 23  --------------   *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                    N       (k1)   (j) (k2)   N'    +-            *
*     ELEMENT:     (j QJ::[ W    * A   ]    ::j QJ)  -+            *
*                                                    ++            *
*                                                    -- B17 (2.3)  *
*                                                                  *
*     SUBROUTINE CALLED: C0T5S,ITJJ3,IXJTIK,IZAS1,RUMTJJ,SIXJ,     *
*                        RMEAJJ,WJ1                                *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      WA=ZERO
      IF(ID(3).EQ.9) THEN
        IF(MAX0(IK(4),ID(4)).LT.3) THEN
          IF(IK(1).LT.300) CALL MES(55)
          IF(ID(1).LT.300) CALL MES(55)
          CALL WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
          RETURN
        ENDIF
      ELSEIF(ID(3).GT.9) THEN
        CALL WAP1JJG(K1,BK2,QM1,QM2,QM3,IK,BK,ID,BD,WA)
        RETURN
      ENDIF
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      KK1=K1*2
      IQ3=QM3*TWO
      KK2=BK2+BK2+TENTH*BK2
      IF(ITJJ3(IK,ID,KK2,BK,BD,IBT,BT,ITP,ITG,IQ3).EQ.0)RETURN
      IQM=TWO*ABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM.GT.IBT(7))GO TO 1
        IF(IXJTIK(KK1,IK(3),KK2,ID(6),IK(6),IBT(6)).EQ.0)GO TO 1
        IBT(1)=IT
        BT(2)=DBLE(IBT(6))/TWO
        BT(1)=DBLE(IBT(7))/TWO
        CALL C0T5S(BD(1),BD(3),QM3,BT(1),BT(3),D1)
        IF(ABS(D1).LT.EPS)GO TO 1
        CALL RMEAJJ(ID(3),IBT(1),IBT(7),IBT(6),ID(1),ID(7),ID(6),S)
        IF(ABS(S).LT.EPS)GO TO 1
        D1=D1*S
        CALL WJ1(IK,BK,IBT,BT,K1,QM1,QM2,W)
        IF(ABS(W).LT.EPS)GO TO 1
        D1=D1*W
        CALL SIXJ(KK1,IK(3),KK2,ID(6),IK(6),IBT(6),0,SI)
        D1=D1*SI/SQRT(DBLE(IBT(7)+1))
        ENQP=ENQP+D1
    1 CONTINUE
      WA=ENQP*SQRT(DBLE(KK2+1))
      IE=KK2+IK(6)+ID(6)+2
      IF(((IE/4)*4).NE.IE)WA=-WA
      RETURN
      END
