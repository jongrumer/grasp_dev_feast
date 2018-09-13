************************************************************************
*                                                                      *
      SUBROUTINE CONSIS (J)
*                                                                      *
*   This routine computes the weighted self-consistency of orbital J   *
*                                                                      *
*   Written by Farid A Parpia, at OXFORD    Last update: 08 Dec 1992   *
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
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /SCF1/UCF(NNNW)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      SCMEA = 0.0D 00
      MTP = MIN (MTP0,MF(J))
      DO 1 I = 1,MTP
         DELTAO = ABS (P(I)-PF(I,J))+ABS (Q(I)-QF(I,J))
         IF (DELTAO .GT. SCMEA) SCMEA = DELTAO
    1 CONTINUE
      SCNSTY(J) = SCMEA*SQRT (UCF(J))
*
      RETURN
      END
