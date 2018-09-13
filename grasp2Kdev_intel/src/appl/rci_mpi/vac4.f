************************************************************************
*                                                                      *
      SUBROUTINE VAC4
*                                                                      *
*   This routine sets up the fourth-order vacuum polarization poten-   *
*   tial using equations (11) and (12) of L Wayne Fullerton and  G A   *
*   Rinker, Jr,  Phys  Rev  A 13 (1976) 1283-1287. The potential  is   *
*   accumulated in array  TC(I), I = 1, ..., N. It is transferred to   *
*   array TA in COMMON block TATB.                                     *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*               [RCI92]: FUNL.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      LOGICAL LDBPR
*
      DIMENSION TC(NNNP)
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF9/CVAC,PI
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NCDIST/ZDIST(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Overall initialization
*
      EPSI = PRECIS*PRECIS
      TWOCV = CVAC+CVAC
*
*   Potential for point nucleus: equation (12)
*
      FACTOR = -Z/(PI*CVAC)**2
*
      TC(1) = 0.0D 00
*
      I = 1
    1 I = I+1
      RI = R(I)
      X = TWOCV*RI
      TCI = (FACTOR/RI)*FUNL (X,1)
      IF (ABS (TCI) .GE. EPSI) THEN
         TC(I) = TCI
         IF (I .LT. N) GOTO 1
      ELSE
         DO 2 K = I,N
            TC(K) = 0.0D 00
    2    CONTINUE
      ENDIF
*
*   Potential for finite nucleus: equation (11)
*
      IF (NPARM .EQ. 2) THEN
*
         FACTOR = -1.0D 00/(PI*CVAC**3)
*
         TC(1) = 0.0D 00
*
         K = 1
    3    K = K+1
*
         RK = R(K)
         XK = TWOCV*RK
         TA(1) = 0.0D 00
*
         DO 4 I = 2,MTP
            XI = TWOCV*R(I)
            XM = ABS (XK-XI)
            XP = XK+XI
            TA(I) = ( FUNL (XM,0)-FUNL (XP,0) )*ZDIST(I)
    4    CONTINUE
*
         CALL QUAD (X)
*
         X = X*FACTOR/RK
*
*   Get out of the loop if the asymptotic region has been reached
*
         IF (ABS (X) .GE. EPSI) THEN
            IF (ABS ((TC(K)-X)/X) .GT. 1.0D-03) THEN
               TC(K) = X
               IF (K .LT. N) GOTO 3
            ENDIF
         ENDIF
*
      ENDIF
*
      IF (LDBPR(8)) THEN
         WRITE (99,300)
         NB2 = N/2
         IF (2*NB2 .EQ. N) THEN
            NROWS = NB2
         ELSE
            NROWS = NB2+1
         ENDIF
         DO 5 II = 1,NROWS
            II1 = II
            II2 = II1+NROWS
            IF (II2 .LE. N) THEN
               WRITE (99,301) R(II1),TB(II1),TC(II1),
     :                        R(II2),TB(II2),TC(II2)
            ELSEIF (II1 .LE. N) THEN
               WRITE (99,301) R(II1),TB(II1),TC(II1)
            ENDIF
    5    CONTINUE
      ENDIF
*
*   Generate total vacuum-polarization potential
*
      DO 6 I = 1,N
         TB(I) = TC(I)+TB(I)
    6 CONTINUE
*
      RETURN
*
  300 FORMAT (///' ++++++++++ VAC4 ++++++++++'
     :         //2(' -------- r -------- ----- VV2 (r) -----'
     :            ,' ----- VV4 (r) -----'))
  301 FORMAT (1P,6(1X,1D19.12))
*
      END
