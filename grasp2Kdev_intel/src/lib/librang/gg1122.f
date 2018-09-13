********************************************************************
*                                                                  *
      SUBROUTINE GG1122(K1,K2,QM1,QM2,QM3,QM4,AA)
*                                                                  *
* ----------------  SECTION METWO    SUBPROGRAM 16  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*     ELEMENTS:                                                    *
*                                                                  *
*       N1       (k1)   N1'         N2       (k2)    N2'        +- *
*     (j  Q J ::W(11)::j  Q'J') * (j  Q J ::W(22):: j  Q'J')    -+ *
*       1  1 1          1  1 1      2  2 2           2  2 2     ++ *
*                                                               -- *
*     SUBROUTINE CALLED: WJ1                                       *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      AA=ZERO
      IF(IK1(3).GT.9) THEN
        IF(IK1(4).GT.2) CALL MES(35)
        IF(ID1(4).GT.2) CALL MES(35)
      ENDIF
      IF(IK2(3).GT.9) THEN
        IF(IK2(4).GT.2) CALL MES(35)
        IF(ID2(4).GT.2) CALL MES(35)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM12=IQMM1+IQMM2
      IF(IK1(4).NE.(ID1(4)+IQMM12))RETURN
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM4=QM4+QM4+TENTH*QM4
      IQMM34=IQMM3+IQMM4
      IF(IK2(4).NE.(ID2(4)+IQMM34))RETURN
      CALL WJ1(IK1,BK1,ID1,BD1,K1,QM1,QM2,A1)
      IF(ABS(A1).LT.EPS)RETURN
      CALL WJ1(IK2,BK2,ID2,BD2,K2,QM3,QM4,W)
      IF(ABS(W).LT.EPS)RETURN
      AA=A1*W
      RETURN
      END
