************************************************************************
*                                                                      *
      SUBROUTINE XPOT (J)
*                                                                      *
*   This subroutine tabulates the exchange terms (the first terms on   *
*   the  right-hand  sides of eqs (14), I P Grant, B J McKenzie, P H   *
*   Norrington, D F Mayers, and N C Pyper, Computer  Phys  Commun 21   *
*   (1980) 211) for orbital J. The exchange terms are stored  in the   *
*   common arrays XP and XQ.                                           *
*                                                                      *
*   Call(s) to: [LIB92]: DRAW, YZK.                                    *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 10 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
Cww      INTEGER PNTNDA,PNTNYA,PNTRDA,PNTRYA
      POINTER (PNTNDA,NDADUMMY), (PNTNYA,NYADUMMY),                     
     :        (PNTRDA,DADUMMY), (PNTRYA,YADUMMY)
      CHARACTER*2 NH
*
      POINTER (PNTRXA,XA(1))
      POINTER (PNTNXA,NXA(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Debug printout: header
*
      IF (LDBPR(27) .OR. LDBPR(28)) WRITE (99,300) NP(J),NH(J)
*
*   Clear for accumulation of sums
*
      DO I = 1, N
         XP(I) = 0.D0
         XQ(I) = 0.D0
      ENDDO
*
*   Add contributions from exchange terms
*
      DO INDEX = 1, NXCOF

         ! Decode information in label
         LABEL = NXA(INDEX)
         K = MOD (LABEL, KEY)
         LABEL = LABEL / KEY
         IOY1  = MOD (LABEL, KEY)
         LABEL = LABEL / KEY
         IOY2  = MOD (LABEL, KEY)
         IORB  = LABEL / KEY
         COEFF = XA(INDEX)

         ! Debug printout: composition
         IF (LDBPR(27)) THEN
            WRITE (99,301) K,COEFF,NP(IOY1),NH(IOY1),NP(IOY2),NH(IOY2),
     :                   NP(IORB),NH(IORB)
         ENDIF

         CALL YZK (K, IOY1, IOY2)
*
*   Accumulate contributions
*
         COEFF = COEFF / C
         !DO I = 1, MF(IORB)
         DO I = 1, n
            CTB = COEFF * TB(I)
            XP(I) = XP(I) + CTB * QF(I,IORB)
            XQ(I) = XQ(I) - CTB * PF(I,IORB)
         ENDDO
      ENDDO
*
*   Debug printout: potential functions
*
      IF (LDBPR(28)) THEN
         WRITE (99,302)
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
               WRITE (99,303) R(II1),XP(II1),XQ(II1),
     :                        R(II2),XP(II2),XQ(II2)
            ELSEIF (II1 .LE. N) THEN
               WRITE (99,303) R(II1),XP(II1),XQ(II1)
            ENDIF
    4    CONTINUE
         CALL DRAW (XP,1.0D 00,XQ,C,N)
      ENDIF
*
      RETURN
*
  300 FORMAT (//' Exchange potential contributions (coefficients will '
     :         ,' be divided by C) for ',1I2,1A2,' orbital :'//)
  301 FORMAT (/25X,'(',1I2,')'
     :        /1X,1PD21.14,'* Y    (',1I2,1A2,',',1I2,1A2,') '
     :                    ,'* P (',1I2,1A2,')')
  302 FORMAT (//31X,'(P)',19X,'(Q)',41X,'(P)',19X,'(Q)'
     :       /2(' --------- r --------- ------ X  (r) -------'
     :         ,' ------ X  (r) -------'))
  303 FORMAT (1P,6(1X,1D21.14))
*
      END
