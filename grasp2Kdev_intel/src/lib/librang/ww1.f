********************************************************************
*                                                                  *
      SUBROUTINE WW1(IK,BK,ID,BD,K2,QM1,QM2,QM3,QM4,WW)
*                                                                  *
*   ---------------  SECTION SQ    SUBPROGRAM 25  --------------   *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                    N      (k2)   (k2) (0)   N'     +-            *
*     ELEMENT      (j QJ::[W   *  W    ]   ::j QJ)   -+            *
*                                                    ++            *
*                                                    -- B17 (2.4)  *
*                                                                  *
*     SUBROUTINE CALLED: ITJJ,IXJTIK,IZAS1,RUMTJJ,WJ1              *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3),IBT(7),BT(3)
      WW=ZERO
      IF(ID(6).NE.IK(6))RETURN
      IF(IZAS1(ID(7),BD(3),IK(7),BK(3)).EQ.0)RETURN
      ENQP=ZERO
      KK2=K2*2
      IQ3=QM3*TWO+QM3*TENTH
      IQ4=QM4*TWO+QM4*TENTH
      IQ=IQ3+IQ4
      IF(ITJJ(IK,ID,0,BK,BD,IBT,BT,KK6,ITP,ITG,IQ).EQ.0)RETURN
      IE1=KK2-IK(6)
      IQM=TWO*ABS(BT(3))+TENTH
      DO 1 IT=ITP,ITG
        CALL RUMTJJ(IT,IBT(3),IBT(7),IBTT,IBT(6))
        IF(IQM.GT.IBT(7))GO TO 1
        IF(IXJTIK(KK2,KK2,0,ID(6),IK(6),IBT(6)).EQ.0)GO TO 1
        IBT(1)=IT
        BT(2)=DBLE(IBT(6))/TWO
        BT(1)=DBLE(IBT(7))/TWO
        CALL WJ1(IK,BK,IBT,BT,K2,QM1,QM2,D1)
        IF(ABS(D1).LT.EPS)GO TO 1
        CALL WJ1(IBT,BT,ID,BD,K2,QM3,QM4,W)
        IF(ABS(W).LT.EPS)GO TO 1
        D1=D1*W
        IE=IE1+IBT(6)
        IF(((IE/4)*4).NE.IE)D1=-D1
        ENQP=ENQP+D1
    1 CONTINUE
      WW=ENQP/SQRT(DBLE(KK2+1)*(IK(6)+1))
      RETURN
      END
