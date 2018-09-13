************************************************************************
*                                                                      *
      FUNCTION SKFUN (K,X)
*                                                                      *
*   Computes the function                                              *
*                                            n  nx                     *
*                             infinity   (-1)  e                       *
*                     S (x) =   Sum      -------                       *
*                               n=1          k                         *
*                                           n                          *
*                                                                      *
*   See, for instance, F. A. Parpia and A. K. Mohanty ``Relativistic   *
*   basis set calculations for atoms with Fermi nuclei'', Phys Rev A   *
*   (1992) in press.                                                   *
*                                                                      *
*   Call(s) to: SKFUN.                                                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
cxhb
cxh place the "magic number" here
      parameter (quasizero=1.d-15)
cxhe

      DNUMER = 1.0D 00
      EN = 0.0D 00
      BASE = -EXP (X)
      SKFUN = 0.0D 00
    1 DNUMER = DNUMER*BASE
      EN = EN+1.0D 00
      DELTA = DNUMER/EN**K
      SKFUN = SKFUN+DELTA
cxhb
cxh      IF (ABS (DELTA/SKFUN) .GT. 1.0D-15 ) GOTO 1
      IF (ABS (DELTA/SKFUN) .GT. quasizero) GOTO 1
cxhe
*
      RETURN
      END
