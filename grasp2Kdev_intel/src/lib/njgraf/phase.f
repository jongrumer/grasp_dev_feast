************************************************************************
*                                                                      *
      SUBROUTINE PHASE (L,JM,NDIM)
*                                                                      *
*   Phase factor arising from non-cyclic permutation of arguments in   *
*   triad L. JM may be either J23 or JDIAG.                            *
*                                                                      *
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
      DIMENSION JM(NDIM,3)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      J7(J7C+1) = JM(L,1)
      J7(J7C+2) = JM(L,2)
      J7C = J7C+3
      J7(J7C) = JM(L,3)
*
      RETURN
      END
