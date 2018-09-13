********************************************************************
*                                                                  *
      SUBROUTINE GG1112(IK1,IK2,BK1,BK2,ID1,ID2,BD1,
     *BD2,K1,QM1,QM2,QM3,QM4,WW)
*                                                                  *
* ----------------  SECTION METWO    SUBPROGRAM 15  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*     ELEMENT:                                                     *
*                                                                  *
*       N1       (k1) (j1)(j2)  N1'       N2     (j2)   N2'     +- *
*     (j Q J ::[W(11)*A(1)]  ::j  Q'J')*(j Q J ::A(2)::j Q'J')  -+ *
*       1 1 1                   1  1 1    2 2 2         2 2 2   ++ *
*                                                               -- *
*     SUBROUTINE CALLED: C0T5S,RMEAJJ,WAP1                         *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK1(7),IK2(7),BK1(3),BK2(3),ID1(7),ID2(7),BD1(3),BD2(3)
      WW=ZERO
      IF(IK1(3).GT.9) THEN
        IF(IK1(4).GT.2) CALL MES(34)
        IF(ID1(4).GT.2) CALL MES(34)
      ENDIF
      IF(IK2(3).GT.9) THEN
        IF(IK2(4).GT.2) CALL MES(34)
        IF(ID2(4).GT.2) CALL MES(34)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM23=IQMM1+IQMM2+IQMM3
      IF(IK1(4).NE.(ID1(4)+IQMM23))RETURN
      IQMM4=QM4+QM4+TENTH*QM4
      IF(IK2(4).NE.(ID2(4)+IQMM4))RETURN
      KK1=K1*2
      CALL C0T5S(BD2(1),BD2(3),QM4,BK2(1),BK2(3),A1)
      IF(ABS(A1).LT.EPS)RETURN
CGG      CALL SJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S)
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),S)
      IF(ABS(S).LT.EPS)RETURN
      BK=HALF*DBLE(IK2(3))
      CALL WAP1(IK1,BK1,ID1,BD1,K1,BK,QM1,QM2,QM3,WA)
      IF(ABS(WA).LT.EPS)RETURN
      WW=-A1*WA*S/SQRT(DBLE(IK2(7)+1))
      RETURN
      END
