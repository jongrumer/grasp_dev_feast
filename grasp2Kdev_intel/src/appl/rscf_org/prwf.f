************************************************************************
*                                                                      *
      SUBROUTINE PRWF (J)
*                                                                      *
*   Makes a (debug) printout of wave functions. There are two modes:   *
*                                                                      *
*      (1) J>0  - Used as a debug option in SOLVE, wavefunctions for   *
*                 orbital J are printed                                *
*      (2) J=0  - A  printout  of the grid and all wave functions is   *
*                 made                                                 *
*                                                                      *
*   Call(s) to: [LIB92]: DRAW.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)

      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      CBZ = C/Z
*
      IF (J .GT. 0) THEN
*
*   Mode (1)
*
         WRITE (99,300) NP(J),NH(J)
         NB2 = MTP0/2
         IF (2*NB2 .EQ. MTP0) THEN
            NROWS = NB2
         ELSE
            NROWS = NB2+1
         ENDIF
         DO 1 II = 1,NROWS
            II1 = II
            II2 = II1+NROWS
            IF (II2 .LE. MTP0) THEN
               WRITE (99,301) R(II1),P(II1),Q(II1),
     :                         R(II2),P(II2),Q(II2)
            ELSEIF (II1 .LE. MTP0) THEN
               WRITE (99,301) R(II1),P(II1),Q(II1)
            ENDIF
    1    CONTINUE
         CALL DRAW (P,1.0D 00,Q,CBZ,MTP0)
*
      ELSE
*
*   Mode (2)
*
         DO 3 K = 1,NW
            WRITE (99,300) NP(K),NH(K)
            MFK = MF(K)
            NB2 = MFK/2
            IF (2*NB2 .EQ. MFK) THEN
               NROWS = NB2
            ELSE
               NROWS = NB2+1
            ENDIF
            DO 2 II = 1,NROWS
               II1 = II
               II2 = II1+NROWS
               IF (II2 .LE. MFK) THEN
                  WRITE (99,301) R(II1),PF(II1,K),QF(II1,K),
     :                            R(II2),PF(II2,K),QF(II2,K)
               ELSEIF (II1 .LE. MFK) THEN
                  WRITE (99,301) R(II1),PF(II1,K),QF(II1,K)
               ENDIF
    2       CONTINUE
            CALL DRAW (PF(1,K),1.0D 00,QF(1,K),CBZ,MF(K))
*
    3    CONTINUE
*
      ENDIF
*
      RETURN
*
  300 FORMAT ('1',1I2,1A2,' orbital:'
     :       /2(' --------- r --------- ------- P (r) -------'
     :         ,' ------- Q (r) -------'))
  301 FORMAT (1P,6(1X,1D21.14))
*
      END
