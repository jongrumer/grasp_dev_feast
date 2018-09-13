************************************************************************
*                                                                      *
      SUBROUTINE SETDM
*                                                                      *
*   Sets dimensions of arrays.                                         *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
*
      JWCC = JWC
      JDELC = JDEL
      J6CC = J6C
      J7CC = J7C
      J8CC = J8C
      J9CC = J9C
*
      RETURN
      END
