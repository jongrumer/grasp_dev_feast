********************************************************************
*                                                                  *
      SUBROUTINE WJ1(IK,BK,ID,BD,K2,QM1,QM2,WJ)
*                                                                  *
*   ---------------  SECTION SQ    SUBPROGRAM 24  --------------   *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                      N       (k2)  N'     +-                     *
*     ELEMENT:       (j  QJ:: W   ::j  QJ)  -+                     *
*                                           ++                     *
*                                           -- S5(1.47),(1.48),    *
*                                                (1.49),(1.50).    *
*                                                                  *
*     SUBROUTINE CALLED: CLE0SM,C1E1SM,RWJJ                        *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),ID(7),BK(3),BD(3)
      WJ=ZERO
      IF(ID(3).EQ.9) THEN
        IF(MAX0(IK(4),ID(4)).LT.3) THEN
          IF(IK(1).LT.300) CALL MES(56)
          IF(ID(1).LT.300) CALL MES(56)
          IQM2=QM2+QM2+QM2*EPS
          IF((ID(4)+IQM2).GT.2) CALL MES(2)
          CALL W1JJG(K2,QM1,QM2,IK,BK,ID,BD,WJ)
          RETURN
        ELSE
          PRINT*, "ERROR in  WJ1"
          STOP
        ENDIF
      ELSEIF(ID(3).GT.9) THEN
        IQM2=QM2+QM2+QM2*EPS
        IF((ID(4)+IQM2).GT.2) CALL MES(2)
        CALL W1JJG(K2,QM1,QM2,IK,BK,ID,BD,WJ)
        RETURN
      ENDIF
      QQ=QM1+QM2
      IF(ABS(QQ).GE.EPS)GO TO 4
      IF(IK(4).NE.ID(4))RETURN
      IF(K2.EQ.0)GO TO 1
      K1=1
      IF(((K2/2)*2).NE.K2)K1=0
      WK1=DBLE(K1)
      CALL CLE0SM(BD(1),BD(3),WK1,BK(1),BK(3),A)
      IF(ABS(A).LT.EPS)RETURN
      CALL RWJJ(IK(3),IK(1),ID(1),K1,K2,W)
      A=A*W
      WJ=A/SQRT(TWO*TWO*BK(1)+TWO)
      IF(QM1.GE.EPS)RETURN
      IF(((K2/2)*2).NE.K2)WJ=-WJ
      RETURN
    1 IF(ID(1).NE.IK(1))RETURN
      IF(QM1.GE.EPS)GO TO 2
      A=DBLE(ID(3)+1-ID(4))
      GO TO 3
    2 A=-DBLE(ID(4))
    3 WJ=A*SQRT(DBLE(IK(6)+1)/DBLE(IK(3)+1))
      RETURN
    4 IF(((K2/2)*2).NE.K2)RETURN
      IQ=QQ+QQ*TENTH
      IF((IK(4)-ID(4)-2*IQ).NE.0)RETURN
      CALL C1E1SM(BD(1),BD(3),QQ,BK(1),BK(3),A)
      IF(ABS(A).LT.EPS)RETURN
      CALL RWJJ(IK(3),IK(1),ID(1),1,K2,W)
      WJ=A*W/SQRT(TWO*BK(1)+ONE)
      RETURN
      END
