************************************************************************
*                                                                      *
      SUBROUTINE BESSEL (IA,IB,IK,IW,K)
*                                                                      *
*   This routine evaluates the functions                               *
*                                                                      *
*                  (2K+1)!!                                            *
*      BESSJ  =    --------  J  (w   r) - 1  =  PHI  (w   r) - 1       *
*                         K   K   ab               K   ab              *
*                  (w   r)                                             *
*                    ab                                                *
*   and                                                                *
*                         K+1                                          *
*                  (w   r)                                             *
*                    ab                                                *
*      BESSN  =  - --------  N  (w   r) - 1  =  PSI  (w   r) - 1       *
*                  (2K-1)!!   K   ab               K   ab              *
*                                                                      *
*                                                                      *
*   where J and N ARE spherical Bessel functions, and                  *
*                                                                      *
*                w   = ABS ( E(IA) - E(IB) )/c                         *
*                 ab                                                   *
*                                                                      *
*   where  E(I)  is the eigenvalue for orbital  I . The writeup (B J   *
*   McKenzie, I P Grant, and P H Norrington, Computer Phys Commun 21   *
*   (1980) 233-246) is incorrect in its description of the output of   *
*   this routine.                                                      *
*                                                                      *
*   The  routine uses equations given in M Abramowitz and I A STegun   *
*   to evaluate the functions. Devices are used to reduce the number   *
*   of actual evaluations of these functions.                          *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      COMMON/BESS1/WIJ(2),BESSJ(2,2,NNNP),BESSN(2,2,NNNP)
     :      /DEBUGR/LDBPR(30)
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STOR/KEEP(2,2)
     :      /WFAC/WFACT
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
      EPSI = SQRT (0.1D 00*ACCY)
*
*   Form unique label symmetric in IA, IB
*
      ICODE = MAX (IA,IB)+KEY*(MIN (IA,IB) + KEY*K )
*
*   Function in position; return
*
      IF (ICODE .EQ. KEEP(IK,IW)) RETURN
*
*   Function not in position; is it available in BESS arrays?
*
      W = WFACT*ABS (E(IA)-E(IB))/C
      WIJ(IW) = W
*
      DO 2 IWKP = 1,2
         DO 2 IKKP = 1,2
            IF ( KEEP(IKKP,IWKP)  .NE.  ICODE) GOTO 2
*
*   Function found move into position
*
            KEEP(IK,IW) = ICODE
            IF (LDBPR(7)) WRITE (99,302) NP(IA),NH(IA),
     :                                   NP(IB),NH(IB),
     :                               K,IKKP,IWKP,IK,IW
            DO 1 L = 1,N
               BESSJ(IK,IW,L) = BESSJ(IKKP,IWKP,L)
               BESSN(IK,IW,L) = BESSN(IKKP,IWKP,L)
    1       CONTINUE
            RETURN
    2 CONTINUE
*
*   Function not found; evaluate it
*
      IF (LDBPR(7)) WRITE (99,303) NP(IA),NH(IA),NP(IB),NH(IB),
     :                              K,IK,IW
*
      KEEP(IK,IW) = ICODE
*
      IF (W .LT. EPSI**2) THEN
*
*   Negligible w
*
         DO 3 I = 1,N
            BESSJ(IK,IW,I) = 0.0D 00
            BESSN(IK,IW,I) = 0.0D 00
    3    CONTINUE
         RETURN
*
      ENDIF
*
      NN = K
*
      BESSJ(IK,IW,1) = 0.0D 00
      BESSN(IK,IW,1) = 0.0D 00
*
*   Use a four-term power series for low w*r
*
      DO 5 J = 2,N
         WA = -0.5D 00*(R(J)*W)**2
         XBESS1 = 1.0D 00
         XBESS2 = 1.0D 00
         S1 = 0.0D 00
         S2 = 0.0D 00
         DO 4 I = 1,4
            XBESS1 = XBESS1*WA/DBLE(I*(2*(NN+I)+1))
            XBESS2 = XBESS2*WA/DBLE(I*(2*(I-NN)-1))
            S1 = S1+XBESS1
            S2 = S2+XBESS2
            IF ((ABS (XBESS1) .LT. ABS (S1)*EPSI) .AND.
     :          (ABS (XBESS2) .LT. ABS (S2)*EPSI)) THEN
               BESSJ(IK,IW,J) = S1
               BESSN(IK,IW,J) = S2
               GOTO 5
            ENDIF
    4    CONTINUE
         JCHAN = J
         GOTO 6
    5 CONTINUE
*
*   If here then calculated whole array using four-term power
*   series.  Hence return
*
      RETURN
*
*   Use sin/cos expansion when power series requires more than
*   four terms terms to converge
*
    6 IF (NN .EQ. 0) THEN
         DFNM = 1.0D 00
         DFN = 1.0D 00
      ELSE
         DFNM = 1.0D 00
         DO 7 I = 3,2*NN-1,2
            DFNM = DFNM*DBLE (I)
    7    CONTINUE
         DFN = DFNM*DBLE (2*NN+1)
      ENDIF
      DFNM = 1.0D 00/DFNM
*
      IREM = MOD (NN,4)
*
      IF (IREM .EQ. 1) THEN
*
*   NN = 1, 5, 9, ...
*
         SSN = -1.0D 00
         SCN =  1.0D 00
         ISWAP = 1
*
      ELSEIF (IREM .EQ. 2) THEN
*
*   N = 2, 6, 10, ....
*
         SSN = -1.0D 00
         SCN = -1.0D 00
         ISWAP = 0
*
      ELSEIF (IREM .EQ. 3) THEN
*
*   NN = 3, 7, 11,...
*
         SSN =  1.0D 00
         SCN = -1.0D 00
         ISWAP = 1
*
      ELSE
*
*   NN = 0, 4, 8,...
*
         SSN =  1.0D 00
         SCN =  1.0D 00
         ISWAP = 0
*
      ENDIF
*
      DO 9 J = JCHAN,N
         WA = W*R(J)
         IF (ISWAP .EQ. 0) THEN
            SN = SSN*SIN (WA)
            CN = SCN*COS (WA)
         ELSE
            SN = SSN*COS (WA)
            CN = SCN*SIN (WA)
         ENDIF
         OBWA = 1.0D 00/WA
         B = OBWA
         S1 = B*SN
         S2 = B*CN
         DO 8 I = 1,NN
            SKEEP = SN
            SN = CN
            CN = -SKEEP
            B =  B*OBWA
     :            *DBLE ((NN+I)*(NN-I+1))
     :            /DBLE (2*I)
            S1 = S1+B*SN
            S2 = S2+B*CN
    8    CONTINUE
         S1 = S1*DFN/(WA**NN)-1.0D 00
         S2 = S2*(WA**(NN+1))*DFNM-1.0D 00
         BESSJ(IK,IW,J) = S1
         BESSN(IK,IW,J) = S2
    9 CONTINUE
      RETURN
*
  303 FORMAT (93X,I2,A2,2X,I2,A2,2X,I2,2X,'New',6X,'(',I2,',',I2,')')
  302 FORMAT (93X,I2,A2,2X,I2,A2,2X,I2,2X,'(',I2,',',I2,')',
     :                                 2X,'(',I2,',',I2,')')
*
      END
