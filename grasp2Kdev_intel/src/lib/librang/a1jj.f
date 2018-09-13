*
*     ---------------------------------------------------------------
*     A 1 J J
*     ---------------------------------------------------------------
*
*
      SUBROUTINE A1JJ(IK,BK,ID,BD,QM1,A)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      DIMENSION IK(7),BK(3),ID(7),BD(3)
      A=ZERO
      IF(QM1.LT.EPS) THEN
        ISUMA=(ID(6)+1)*ID(4)
        AB=DBLE(ISUMA)
        A=DSQRT(AB)
        IFAZ=ID(6)+ID(3)-IK(6)+ID(4)*2
        IF((IFAZ/4)*4.NE.IFAZ)A=-A
      ELSE
        ISUMA=(IK(6)+1)*IK(4)
        AB=DBLE(ISUMA)
        A=DSQRT(AB)
        IF((IK(4)/2)*2.NE.IK(4))A=-A
      ENDIF
      RETURN
      END
