************************************************************************
*                                                                      *
      SUBROUTINE DAMPCK (IPR,J,ED1,ED2)
*                                                                      *
*   This subroutine determines the damping factor appropriate to the   *
*   present  orbital. The algorithm is taken from C Froese Fischer's   *
*   program MCHF.                                                      *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
*   Modified by C. Froese Fischer           Last update: 07 Apr 2009   *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PCDAMP
      POINTER (PCDAMP,DAMPDUMMY)
      LOGICAL ADAPTV
*
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /ORB1/E(NNNW),GAMA(NNNW),PED(NNNW)
*
*   The damping is adaptive (i.e., can be modified by this SUBROUTINE)
*   if and only if ODAMP(J) .GE. 0.0
*
      ADAPTV = ODAMP(J) .GE. 0.0D 00
*
      ED2 = (ED2-E(J))/ED2
*     Print *, 'ED1,ED2,ED1*ED2', ed1,ed2,ed1*ed2
      IF (ADAPTV) THEN
         IF (ED1*ED2 .LT. -1.0d-4) THEN
*           We have a significant oscillation, damp
            ODAMP(J) = 0.10D 00+0.90D 00*ODAMP(J)
*           Print *, 'Increasing ODAMP', odamp(j)
         ELSE 
               ODAMP(J) = 0.50D 00*ODAMP(J)
         ENDIF
      ENDIF
*
*     Save the relative difference for the next update of this orbital
      PED(J) = ED2
*
      RETURN
      END
