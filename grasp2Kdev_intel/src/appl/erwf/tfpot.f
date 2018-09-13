************************************************************************
*                                                                      *
      SUBROUTINE TFPOT
*                                                                      *
*   Calculation of the universal Thomas-Fermi potential.               *
*                                                                      *
*   Call(s) to: [LIB92]: DRAW.                                         *
*                                                                      *
*                                         Last revision: 09 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL LDBPR
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
*
      THIRD = 1.0D 00/3.0D 00
*
      WA = Z-DBLE (NELEC-1)
      WB = MAX (Z-WA,0.0D 00)
      WB = WB**THIRD/0.8853D 00
      DO 1 I = 1,N
*
*   Rational function approximation to the universal Thomas-Fermi
*   function
*
         WC = SQRT (WB*R(I))
         WD = WC*(0.60112D0*WC+1.81061D0)+1.0D 00
         WE = WC*(WC*(WC*(WC*(0.04793D0*WC+0.21465D0)+
     :      0.77112D0)+1.39515D0)+1.81061D0)
     :      +1.0D 00
         WF = WD/WE
         YP(I) = (ZZ(I)-WA)*WF*WF+WA
    1 CONTINUE
*
*   Debug printout
*
      IF (LDBPR(26)) THEN
         WRITE (99,300)
         NB3 = N/3
         IF (3*NB3 .EQ. N) THEN
            NROWS = NB3
         ELSE
            NROWS = NB3+1
         ENDIF
         DO 3 II = 1,NROWS
            II1 = II
            II2 = II1+NROWS
            II3 = II2+NROWS
            IF (II3 .LE. N) THEN
               WRITE (99,301) R(II1),YP(II1),R(II2),YP(II2),
     :                         R(II3),YP(II3)
            ELSEIF (II2 .LE. N) THEN
               WRITE (99,301) R(II1),YP(II1),R(II2),YP(II2)
            ELSE
               WRITE (99,301) R(II1),YP(II1)
            ENDIF
    3    CONTINUE
         CALL DRAW (YP,1.0D 00,YP,1.0D 00,N)
      ENDIF
*
      RETURN
*
  300 FORMAT (///' Thomas-Fermi potential'
     :       //3(' --------- r --------- ------ -r*V(r) ------'))
  301 FORMAT (1P,6(1X,1D21.14))
*
      END
