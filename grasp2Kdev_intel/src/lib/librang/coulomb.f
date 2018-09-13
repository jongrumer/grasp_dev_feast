********************************************************************
*                                                                  *
      SUBROUTINE COULOM(J1,J2,J3,J4,L1,L2,L3,L4,K,AA)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 01  ------------   *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF COULOMB INTERACTIONS BETWEEN THE ELECTRONS                *
*                                                                  *
*                          k   k+1  (k) (k)                        *
*     (n l j T  n l j T ::r  / r  ( C   C )::n l j T  n l j T )    *
*       1 1 1 1  2 2 2 2   <    >             3 3 3 3  4 4 4 4     *
*                                                                  *
*     SUBROUTINE CALLED:  CRE                                      *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      AA=ZERO
      IF(ITTK(L1,L3,K).EQ.0)RETURN
      IF(ITTK(L2,L4,K).EQ.0)RETURN
      I=(2*K+1)/2
      AA=CRE (J1,I,J3)
      IF(ABS(AA).LT.EPS)RETURN
      AA=AA*CRE (J2,I,J4)
      IF(ABS(AA).LT.EPS)RETURN
      IFAZ=L3-2*K-L1+L4-L2
      IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
      RETURN
      END
