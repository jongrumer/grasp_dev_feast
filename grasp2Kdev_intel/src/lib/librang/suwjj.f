*
*     -------------------------------------------------------------
*     S U W J J
*     -------------------------------------------------------------
*
*                 (k1 k2)
*     ( j QJ ::: W      ::: j QJ )
*
*
      SUBROUTINE SUWJJ(K1,K2,LL,J1,J2,SUW)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/RIBOJJ/IMPTJJ(63),IMGTJJ(63),IMPNJJ(63),IMGNJJ(63)
      SUW=ZERO
      IF(IMPTJJ(J1).NE.IMPTJJ(J2)) RETURN
      S=ZERO
      CALL RUMTJJ(J1,LL,L2Q1,L2V1,L2J1)
      KK1=K1*2
      KK2=K2*2
      IP=IMPNJJ(J1)
      IG=IMGNJJ(J1)
      CALL RUMTJJ(J2,LL,L2Q2,L2V2,L2J2)
      DO 1 I=IP,IG
        CALL RUMTJJ(I,LL,L2QI,L2VI,L2JI)
        IF(IXJTIK(LL,LL,KK2,L2J2,L2J1,L2JI).NE.0) THEN
	    IF(IXJTIK(1,1,KK1,L2Q2,L2Q1,L2QI).NE.0) THEN
            CALL RMEAJJ(LL,J1,L2Q1,L2J1,I,L2QI,L2JI,COEF1)
            IF(DABS(COEF1).GT.EPS) THEN
              CALL RMEAJJ(LL,I,L2QI,L2JI,J2,L2Q2,L2J2,COEF2)
              IF(DABS(COEF2).GT.EPS) THEN
                CALL SIXJ(LL,LL,KK2,L2J2,L2J1,L2JI,0,SI1)
                CALL SIXJ(1,1,KK1,L2Q2,L2Q1,L2QI,0,SI2)
                S=S+SI1*SI2*COEF1*COEF2
              ENDIF
            ENDIF
          ENDIF
        ENDIF
    1 CONTINUE
      SUW=S*DSQRT(DBLE((KK1+1)*(KK2+1)))
      IF(MOD(L2Q1+L2J1+L2Q2+L2J2+KK1+KK2,4).NE.0)SUW=-SUW
      RETURN
      END
