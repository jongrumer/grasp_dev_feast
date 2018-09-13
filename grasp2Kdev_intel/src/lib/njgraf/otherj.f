************************************************************************
*                                                                      *
      SUBROUTINE OTHERJ (LIN,J,LO,LCO,K)
*                                                                      *
*   Gives the other triad where a given J occurs and its position.     *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INTEGER ARROW
      LOGICAL TABS
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
*
      LO = LINE(J,1)
      IF ((LO .EQ. LIN) .OR. (TABS(LO))) THEN
         K = 1
         LO = LINE(J,2)
         LCO = LCOL(J,2)
      ELSE
         K = 2
         LCO = LCOL(J,1)
      ENDIF
*
      RETURN
      END
