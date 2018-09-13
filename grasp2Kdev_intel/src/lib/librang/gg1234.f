********************************************************************
*                                                                  *
      SUBROUTINE GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,
     *ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,WW)
*                                                                  *
* ----------------  SECTION METWO    SUBPROGRAM 19  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                      N1     (j1)  N1'      N2    (j2)   N2'      *
*     ELEMENTS:      (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  *
*                      1 1 1        1 1 1    2 2 2        2 2 2    *
*                                                                  *
*        N3     (j3)  N3'      N4     (j4)  N4'                 +- *
*     *(j Q J::A(3)::j Q'J')*(j Q J::A(4)::j Q'J')              -+ *
*        3 3 3        3 3 3    4 4 4        4 4 4               ++ *
*                                                               -- *
*     SUBROUTINE CALLED: C0T5S,RMEAJJ                              *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK1(7),IK2(7),IK3(7),IK4(7),BK1(3),BK2(3),BK3(3),
     *BK4(3),ID1(7),ID2(7),ID3(7),ID4(7),BD1(3),BD2(3),
     *BD3(3),BD4(3)
      WW=ZERO
      IF(IK1(3).GT.9) THEN
        IF(IK1(4).GT.2) CALL MES(30)
        IF(ID1(4).GT.2) CALL MES(30)
      ENDIF
      IF(IK2(3).GT.9) THEN
        IF(IK2(4).GT.2) CALL MES(30)
        IF(ID2(4).GT.2) CALL MES(30)
      ENDIF
      IF(IK3(3).GT.9) THEN
        IF(IK3(4).GT.2) CALL MES(30)
        IF(ID3(4).GT.2) CALL MES(30)
      ENDIF
      IF(IK4(3).GT.9) THEN
        IF(IK4(4).GT.2) CALL MES(30)
        IF(ID4(4).GT.2) CALL MES(30)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4).NE.(ID1(4)+IQMM1))RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IF(IK2(4).NE.(ID2(4)+IQMM2))RETURN
      IQMM3=QM3+QM3+TENTH*QM3
      IF(IK3(4).NE.(ID3(4)+IQMM3))RETURN
      IQMM4=QM4+QM4+TENTH*QM4
      IF(IK4(4).NE.(ID4(4)+IQMM4))RETURN
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
      IF(ABS(A1).LT.EPS)RETURN
      CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
      IF(ABS(C1).LT.EPS)RETURN
      A1=A1*C1
      CALL C0T5S(BD3(1),BD3(3),QM3,BK3(1),BK3(3),C1)
      IF(ABS(C1).LT.EPS)RETURN
      A1=A1*C1
      CALL C0T5S(BD4(1),BD4(3),QM4,BK4(1),BK4(3),C1)
      IF(ABS(C1).LT.EPS)RETURN
      A1=A1*C1
CGG      CALL SJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      IF(ABS(S).LT.EPS)RETURN
CGG      CALL SJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
      IF(ABS(C).LT.EPS)RETURN
CGG      CALL SJJ(IK3(3),IK3(1),IK3(7),IK3(6),ID3(1),ID3(7),ID3(6),V)
      CALL RMEAJJ(IK3(3),IK3(1),IK3(7),IK3(6),ID3(1),ID3(7),ID3(6),V)
      IF(ABS(V).LT.EPS)RETURN
CGG      CALL SJJ(IK4(3),IK4(1),IK4(7),IK4(6),ID4(1),ID4(7),ID4(6),Z)
      CALL RMEAJJ(IK4(3),IK4(1),IK4(7),IK4(6),ID4(1),ID4(7),ID4(6),Z)
      IF(ABS(Z).LT.EPS)RETURN
      WW=A1*S*C*V*Z/SQRT(DBLE((IK1(7)+1)*(IK2(7)+1)*
     *(IK3(7)+1)*(IK4(7)+1)))
      RETURN
      END
