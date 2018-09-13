************************************************************************
*                                                                      *
      SUBROUTINE SETXZ (J)
*                                                                      *
*   This subprogram sets the inhomogeneous terms to zero.              *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
*
      DO 1 I = 1,N
         XU(I) = 0.0D 00
         XV(I) = 0.0D 00
    1 CONTINUE
*
      RETURN
      END
