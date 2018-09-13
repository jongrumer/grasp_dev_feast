************************************************************************
*                                                                      *
      SUBROUTINE LAGCON (J, nprocs)
*                                                                      *
*   This  routine  includes  the Lagrange multiplier contribution in   *
*   the 'exchange' term.                                               *
*
*   Parameter nprocs is added so that it works for mpi and serial 
*   programs.
*                                                                      *
*                                           Last update: 08 Dec 1992   *
*   Modified by Xinghong He                 Last update: 17 Aug 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      POINTER (PNTECV,ECV(1))
      POINTER (PNIECC,IECC(1))
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /LAGR/PNTECV,PNIECC,NEC
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF1/UCF(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
      IF (NEC .EQ. 0) RETURN
*
Cww      EPS = ACCY*0.1D 00
      EPS = ACCY*0.01D 00
*
*   Add contributions from off-diagonal parameters to exchange
*
      WB = 1.0D 00/(UCF(J)*C) / nprocs
      DO 2 K = 1,NEC
*
*   Decode index
*
         IECCK = IECC(K)
         L1 = IECCK/KEY
         L2 = IECCK-KEY*L1
*
         IF (J .EQ. L1) THEN
            M = L2
         ELSEIF (J .EQ. L2) THEN
            M = L1
         ELSE
            GOTO 2
         ENDIF
*
         WA = ECV(K)*WB
         IF (ABS (WA) .LT. EPS) GOTO 2
*
*   ADD CONTRIBUTIONS TO EXCHANGE TERMS
*
         MFM = MF(M)
         DO 1 I = 1,MFM
            WARI = WA*R(I)
            XP(I) = XP(I)+WARI*QF(I,M)
            XQ(I) = XQ(I)-WARI*PF(I,M)
    1    CONTINUE
*
    2 CONTINUE
*
      RETURN
      END
