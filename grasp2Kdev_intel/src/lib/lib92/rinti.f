************************************************************************
*                                                                      *
      FUNCTION RINTI (J,K,MODE)
*                                                                      *
*   The value of this  function is the one-electron integral I (J,K)   *
*   for  orbitals  J, K. The analytical expression for this quantity   *
*   is given as  eq (9) in  I P Grant, B J McKenzie, P H Norrington,   *
*   D F Mayers, and N C Pyper,  Computer  Phys Commun 21 (1980) 211 .  *
*                                                                      *
*   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 06 Oct 1992   *
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
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Stop if orbitals J and K have different kappa values
*
      IF (NAK(J) .NE. NAK(K)) THEN
         WRITE (*,300) NP(J),NH(J),NP(K),NH(K)
         STOP
      ENDIF
*
      MTP = MAX (MF(J),MF(K))
*
*   Kinetic energy contribution
*
*   Piece involving derivatives
*
      CALL DPBDT (K)
      TA(1) = 0.D0
      DO 1 I = 2,MTP
         TA(I) = QF(I,J)*TA(I)-PF(I,J)*TB(I)
    1 CONTINUE
      CALL QUAD (PIECE1)
      PIECE1 = C*PIECE1/H
*
*   Pieces not involving derivatives
*
      TA(1) = 0.D0
      DO 2 I = 2,MTP
         TA(I) = RP(I)*QF(I,J)*QF(I,K)
    2 CONTINUE
      CALL QUAD (PIECE2)
      PIECE2 = -2.D0*C*C*PIECE2
*
      TA(1) = 0.D0
      DO 3 I = 2,MTP
         TA(I) = RPOR(I)*(PF(I,J)*QF(I,K)+QF(I,J)*PF(I,K))
    3 CONTINUE
      CALL QUAD (PIECE3)
      PIECE3 = PIECE3*C*DBLE (NAK(K))
*
*   Contribution from nuclear potential only if MODE is 0
*
      IF (MODE .EQ. 0) THEN
*
         TA(1) = 0.D0
         DO 4 I = 2,MTP
            TA(I) = RPOR(I)*ZZ(I)
     :                     *(PF(I,J)*PF(I,K)+QF(I,J)*QF(I,K))
    4    CONTINUE
         CALL QUAD (PIECE4)
         PIECE4 = -PIECE4
*
      ELSE
         PIECE4 = 0.D0
      ENDIF
*
      RINTI = PIECE1+PIECE2+PIECE3+PIECE4
*
*   Debug printout
*
      IF ((MODE .EQ. 0) .AND. LDBPR(4))
     :   WRITE (99,301) NP(J),NH(J),NP(K),NH(K),RINTI
      IF ((MODE .NE. 0) .AND. LDBPR(5))
     :   WRITE (99,302) NP(J),NH(J),NP(K),NH(K),RINTI
*
      RETURN
*
  300 FORMAT ('RINTI: Attempt to calculate I(',1I2,1A2,',',1I2,1A2,')')
  301 FORMAT (/' I (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12/)
  302 FORMAT (/' K (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12/)
*
      END
