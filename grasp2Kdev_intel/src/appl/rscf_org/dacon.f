************************************************************************
*                                                                      *
      SUBROUTINE DACON
*                                                                      *
*   This  routine  includes  the  contribution from the off-diagonal   *
*   I(a,b) integrals in the 'exchange' term.                           *
*                                                                      *
*   Call(s) to: [LIB92]: DPBDT.                                        *
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

      POINTER (PNTNXA,NXADUMMY)
      POINTER (PNTNYA,NYADUMMY)
      POINTER (PNTRXA,RXADUMMY)
      POINTER (PNTRYA,RYADUMMY)
*
      POINTER (PNTRDA,DA(1))
      POINTER (PNTNDA,NDA(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPOT/ZZ(NNNP),NNUC
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF1/UCF(NNNW)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      TWOC = C + C
*
      DO 2 K = 1, NDCOF
*
         IORB = NDA(K)
         CALL DPBDT (IORB)
         MFI = MF(IORB)
         COEFF = DA(K)
         FK = DBLE (NAK(IORB))
*
         DO 1 I = 2, MFI
            RPORII = 1.0D0/(H*RPOR(I))
            PFI = PF(I,IORB)
            QFI = QF(I,IORB)
            ZBCI = ZZ(I)/C
            XP(I) = XP(I)+COEFF*( TA(I)*RPORII+FK*PFI
     :                           -(TWOC*R(I)+ZBCI)*QFI )
            XQ(I) = XQ(I)+COEFF*( TB(I)*RPORII-FK*QFI
     :                                      +ZBCI *PFI )
    1    CONTINUE
*
    2 CONTINUE
*
      RETURN
      END
