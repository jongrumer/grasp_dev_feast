********************************************************************
*                                                                  *
      SUBROUTINE GG12(IK1,IK2,BK1,BK2,ID1,ID2,BD1,BD2,QM1,QM2,WW)
*                                                                  *
* ----------------  SECTION METWO    SUBPROGRAM 18  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF FOLLOWING MATRIX       *
*                                                                  *
*                   N1     (j1)  N1'      N2     (j2)  N2'      +- *
*     ELEMENTS:   (j Q J::A(1)::j Q'J')*(j Q J::A(2)::j Q'J')*  -+ *
*                   1 1 1        1 1 1    2 2 2        2 2 2    ++ *
*                                                               -- *
*     SUBROUTINE CALLED: C0T5S,RMEAJJ,WJ1                          *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK1(7),IK2(7),BK1(3),BK2(3),ID1(7),ID2(7),BD1(3),BD2(3)
      WW=ZERO
      IF(IK1(3).GT.9) THEN
        IF(IK1(4).GT.2) CALL MES(32)
        IF(ID1(4).GT.2) CALL MES(32)
      ENDIF
      IF(IK2(3).GT.9) THEN
        IF(IK2(4).GT.2) CALL MES(32)
        IF(ID2(4).GT.2) CALL MES(32)
      ENDIF
      IQMM1=QM1+QM1+TENTH*QM1
      IF(IK1(4).NE.(ID1(4)+IQMM1))RETURN
      IQMM2=QM2+QM2+TENTH*QM2
      IF(IK2(4).NE.(ID2(4)+IQMM2))RETURN
      CALL C0T5S(BD1(1),BD1(3),QM1,BK1(1),BK1(3),A1)
      IF(ABS(A1).LT.EPS)RETURN
      CALL C0T5S(BD2(1),BD2(3),QM2,BK2(1),BK2(3),C1)
      IF(ABS(C1).LT.EPS)RETURN
      A1=A1*C1
      CALL RMEAJJ(IK1(3),IK1(1),IK1(7),IK1(6),ID1(1),ID1(7),ID1(6),S)
      IF(ABS(S).LT.EPS)RETURN
      CALL RMEAJJ(IK2(3),IK2(1),IK2(7),IK2(6),ID2(1),ID2(7),ID2(6),C)
      IF(ABS(C).LT.EPS)RETURN
      WW=A1*S*C/SQRT(DBLE((IK1(7)+1)*(IK2(7)+1)))
      RETURN
      END
