************************************************************************
*                                                                      *
      SUBROUTINE SETXV (J)
*                                                                      *
*   This  subprogram  sets up the inhomogeneous terms for the varia-   *
*   tion equations.                                                    *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
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
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      HHC = 0.5D 00*H/C
*
*   Set up arrays TF and TG
*
      DO 1 I = 1,N
         XU(I) = -QF(I,J)*HHC*RP(I)
         XV(I) =  PF(I,J)*HHC*RP(I)
    1 CONTINUE
*
      DO 2 I = 2,N
         XU(I-1) = XU(I)+XU(I-1)
         XV(I-1) = XV(I)+XV(I-1)
    2 CONTINUE
*
      RETURN
*
      END
