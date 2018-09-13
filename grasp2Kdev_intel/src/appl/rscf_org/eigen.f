************************************************************************
*                                                                      *
      FUNCTION EIGEN (J)
*                                                                      *
*   This function computes an estimate of the energy of orbital J .    *
*                                                                      *
*   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Initialization
*
      MTP = MF(J)
*
*   Exchange term
*
      TA(1) = 0.0D 00
      DO 1 I = 2,MTP
         TA(I) = (PF(I,J)*XQ(I)-QF(I,J)*XP(I))*RPOR(I)
    1 CONTINUE
      CALL QUAD (PIECE1)
      PIECE1 = C*PIECE1
*
*   Direct term
*
      TA(1) = 0.0D 00
      DO 2 I = 2,MTP
         TA(I) = (PF(I,J)**2+QF(I,J)**2)*YP(I)*RPOR(I)
    2 CONTINUE
      CALL QUAD (PIECE2)
*
*   Kinetic energy terms
*
      TA(1) = 0.0D 00
      DO 3 I = 2,MTP
         TA(I) = (QF(I,J)**2)*RP(I)
    3 CONTINUE
      CALL QUAD (PIECE3)
      PIECE3 = 2.0D 00*C*C*PIECE3
*
      TA(1) = 0.0D 00
      DO 4 I = 2,MTP
         TA(I) = (PF(I,J)*QF(I,J))*RPOR(I)
    4 CONTINUE
      CALL QUAD (PIECE4)
      PIECE4 = -2.0D 00*DBLE (NAK(J))*C*PIECE4
*
      CALL DPBDT (J)
      TA(1) = 0.0D 00
      DO 5 I = 2,MTP
         TA(I) = PF(I,J)*TB(I)-QF(I,J)*TA(I)
    5 CONTINUE
      CALL QUAD (PIECE5)
      PIECE5 = C*PIECE5/H
*
*   Assembly
*
      EIGEN = PIECE1+PIECE2+PIECE3+PIECE4+PIECE5
*
      RETURN
*
      END
