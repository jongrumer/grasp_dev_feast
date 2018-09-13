************************************************************************
*                                                                      *
      SUBROUTINE RADGRD
*                                                                      *
*   This routine sets up the radial grid  R  and the associated arr-   *
*   ays  RP  and  RPOR  in the COMMON block  /GRID/. Different grids   *
*   are generated depending on whether or not  HP  is  zero .          *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      LOGICAL LDBPR
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
*
*   RPOR(1) is never used in the program: it is arbitrarily
*   set to zero
*
      NP10 = N+10
      R(1) = 0.0D 00
      RPOR(1) = 0.0D 00
*
*   Now set up the grids
*
      IF (HP .EQ. 0.0D 00) THEN
!        default comes here
*
*   Exponential grid if HP is zero
*
*   Initializations
*
         RP(1) = RNT
         EPH = EXP (H)
         ETT = 1.0D 00
*
*   Set up the arrays R, RP, RPOR
*
         DO 1 I = 2,NP10
            ETT = EPH*ETT
            ETTM1 = ETT-1.0D 00
            R(I) = RNT*ETTM1
            RP(I) = RNT*ETT
            RPOR(I) = ETT/ETTM1
    1    CONTINUE
*
      ELSE
*
*   Asymptotically-linear exponential grid otherwise:
*
*   Initializations
*
         EPSLON = 1.0D 03*PRECIS
         A = H/HP
         RP(1) = RNT/(A*RNT+1.0D 00)
         RLAST = 0.0D 00
         REST = 0.0D 00
*
*   Set up the arrays R, RP, RPOR
*
         DO 3 I = 2,NP10
*
            T = H*DBLE (I-1)
*
*   Solve the implicit equation for R using the Newton-Raphson
*   method
*
    2       RESTS = REST+RNT
            FOFR = LOG (RESTS/RNT)+A*REST-T
            FPRI = RESTS/(A*RESTS+1.0D 00)
            DELR = -FOFR*FPRI
            REST = RLAST+DELR
*
            IF (ABS (DELR/REST) .LT. EPSLON) THEN
               R(I) = REST
               RESTS = REST+RNT
               FPRI = RESTS/(A*RESTS+1.0D 00)
               RP(I) = FPRI
               RPOR(I) = FPRI/REST
            ELSE
               RLAST = REST
               GOTO 2
            ENDIF
*
    3    CONTINUE
*
      ENDIF
*
*   Debug printout
*
      IF (LDBPR(1)) THEN
         WRITE (99,300)
         NB2 = N/2
         IF (2*NB2 .EQ. N) THEN
            NROWS = NB2
         ELSE
            NROWS = NB2+1
         ENDIF
         DO 4 II = 1,NROWS
            II1 = II
            II2 = II1+NROWS
            IF (II2 .LE. N) THEN
               WRITE (99,301) R(II1),RP(II1),RPOR(II1),
     :                         R(II2),RP(II2),RPOR(II2)
            ELSEIF (II1 .LE. N) THEN
               WRITE (99,301) R(II1),RP(II1),RPOR(II1)
            ENDIF
    4    CONTINUE
      ENDIF
*
      RETURN
*
  300 FORMAT (/'From SUBROUTINE RADGRD:'
     :        /2(' -------- r -------- -------- r'' -------'
     :          ,' ------- r''/r ------'))
  301 FORMAT (1P,6(1X,1D19.12))
*
      END
