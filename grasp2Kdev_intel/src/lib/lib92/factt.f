************************************************************************
*                                                                      *
      SUBROUTINE FACTT
*                                                                      *
*   Calculates the logs  of factorials required by the Racah coeffi-   *
*   cient routine DRACAH.                                              *
*                                                                      *
*   Written by N S Scott                    Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      PARAMETER (MFACT = 500)
*
      COMMON/FACTS/GAM(MFACT)
*
      GAM(1) = 1.0D 00
      GAM(2) = 1.0D 00
      X = 2.0D 00
*
      DO 10 I = 3,30
         GAM(I) = GAM(I-1)*X
         X = X+1.0D 00
   10 CONTINUE
*
      DO 20 I = 1,30
         GAM(I) = LOG(GAM(I))
  20  CONTINUE
*
      X = 3.0D 01
*
      DO 30 I = 31,MFACT
         GAM(I) = GAM(I-1)+LOG(X)
         X = X+1.0D 00
  30  CONTINUE
*
      RETURN
      END
