************************************************************************
*                                                                      *
      SUBROUTINE ES (F,S2F,S3f)
*                                                                      *
*   Evaluate the sum of the series                                     *
*                                                                      *
*                       infinity      n              k                 *
*              S  (F) =   Sum     (-1)  exp (n*F) / n                  *
*               k        n = 0                                         *
*                                                                      *
*   for k = 2, 3 to machine precision. This is a utility subroutine,   *
*   called by SUBROUTINEs NUCPOT and NCHARG.                           *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 28 Sep 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      N = 0
      S2F = 0.0D 00
      S3F = 0.0D 00
      FASE = 1.0D 00
    1 N = N+1
         EN = DBLE (N)
         OBN = 1.0D 00/EN
         FASE = -FASE
         ENF = EXP (EN*F)
         TERM2 = FASE*ENF*OBN*OBN
         TERM3 = TERM2*OBN
         S2LAST = S2F
         S2F = S2F+TERM2
         S3F = S3F+TERM3
         IF (ABS (S2F) .NE. ABS (S2LAST)) GOTO 1
*
      RETURN
      END
