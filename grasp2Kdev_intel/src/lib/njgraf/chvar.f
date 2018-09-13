************************************************************************
*                                                                      *
      SUBROUTINE CHVAR (JP,NBC,KBC,JT,JINV,NSUM)
*                                                                      *
*   Change  the  order  of  summation variable to be able to perform   *
*   separately the summations in GENSUM.                               *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL JT
*
      DIMENSION JINV(NSUM),JP(NBC),JT(NSUM)
*
      KB = KBC+1
*
      IF (KB .LE. NBC) THEN
         DO 1 I = KB,NBC
            JK = JP(I)
            IF (JT(JK)) THEN
               KBC = KBC+1
               JP(I) = JP(KBC)
               JP(KBC) = JINV(JK)
            ENDIF
    1    CONTINUE
      ENDIF
*
      RETURN
      END
