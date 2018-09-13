************************************************************************
*                                                                      *
      FUNCTION CRE (KAP1,K,KAP2)
*                                                                      *
*   Computes the relativistic reduced matrix element                   *
*                                                                      *
*                         (j1 || C(K) || j2),                          *
*                                                                      *
*   Eq. (5.15) of I P Grant, Advances in Physics 19 (1970) 762. KAP1,  *
*   KAP2 are the kappa values corresponding to j1, j2.  The triangle   *
*   conditions are tested by the routine CLRX.                         *
*                                                                      *
*   Call(s) to: [LIB92] CLRX.                                          *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      K1 = ABS (KAP1)
      DK1K2 = DBLE (4*K1*IABS (KAP2))
      CRE = SQRT (DK1K2)*CLRX (KAP1,K,KAP2)
      IF (MOD (K1,2) .EQ. 1) CRE  = -CRE
*
      RETURN
      END
