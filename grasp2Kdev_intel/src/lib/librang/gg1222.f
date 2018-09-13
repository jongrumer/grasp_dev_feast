********************************************************************
*                                                                  *
      SUBROUTINE GG1222(IK1,IK2,BK1,BK2,ID1,ID2,BD1,
     *BD2,K1,QM1,QM2,QM3,QM4,WW)
*                                                                  *
* ----------------  SECTION METWO    SUBPROGRAM 17  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                                          N1      (j1)   N1'      *
*     ELEMENTS:                          (j  Q J ::A(1)::j Q'J')*  *
*                                          1  1 1         1 1 1    *
*        N2        (j2)    (k1) (j1)  N2'                       +- *
*     *(j  Q J ::[ A(2) * W(22) ]  ::j  Q'J')                   -+ *
*        2  2 2                       2  2 2                    ++ *
*                                                               -- *
*     SUBROUTINE CALLED: C0T5S,RMEAJJ,AWP1                         *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK1(7),IK2(7),BK1(3),BK2(3),ID1(7),ID2(7),
     *BD1(3),BD2(3)
      WW=ZERO
      IF(IK1(3).GT.9) THEN
        IF(IK1(4).GT.2) CALL MES(33)
        IF(ID1(4).GT.2) CALL MES(33)
      ENDIF
      IF(IK2(3).GT.9) THEN
        IF(IK2(4).GT.2) CALL MES(33)
        IF(ID2(4).GT.2) CALL MES(33)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4).NE.(ID1(4)+IQMM1))RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IQMM3=QM3+QM3+TENTH*QM3
      IQMM4=QM4+QM4+TENTH*QM4
      IQMM34=IQMM2+IQMM3+IQMM4
      IF(IK2(4).NE.(ID2(4)+IQMM34))RETURN
      KK1=K1*2
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
      IF(ABS(A1).LT.EPS)RETURN
CGG      CALL SJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      IF(ABS(S).LT.EPS)RETURN
      BK=HALF*DBLE(IK1(3))
      CALL AWP1(IK2,BK2,ID2,BD2,K1,BK,QM2,QM3,QM4,AW)
      IF(ABS(AW).LT.EPS)RETURN
      WW=-A1*AW*S/SQRT(DBLE(IK1(7)+1))
      RETURN
      END
