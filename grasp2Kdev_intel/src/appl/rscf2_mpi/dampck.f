************************************************************************
*                                                                      *
      SUBROUTINE DAMPCK (IPR,J,ED1,ED2)
*                                                                      *
*   This subroutine determines the damping factor appropriate to the   *
*   present  orbital. The algorithm is taken from C Froese Fischer's   *
*   program MCHF.                                                      *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PCDAMP
      POINTER (PCDAMP,DAMPDUMMY)
      LOGICAL ADAPTV
*
      COMMON/DAMP/ODAMP(NNNW),PCDAMP
     :      /ORB1/E(NNNW),GAMA(NNNW)
*
*   The damping is adaptive (i.e., can be modified by this SUBROUTINE)
*   if and only if ODAMP(J) .GE. 0.0
*
      ADAPTV = ODAMP(J) .GE. 0.0D 00
*
      IF (IPR .NE. J) THEN
         IF (ADAPTV) ODAMP(J) = 0.75D 00*ODAMP(J)
      ELSE
         ED2 = ED2-E(J)
         IF (ADAPTV) THEN
            IF (ED1*ED2 .GT. 0.0D 00) THEN
               ODAMP(J) = 0.75D 00*ODAMP(J)
            ELSE
               ODAMP(J) = 0.25D 00+0.75D 00*ODAMP(J)
            ENDIF
         ENDIF
      ENDIF
*
      IF (IPR .EQ. J) THEN
         ED1 = ED2
      ELSE
         ED1 = ED2-E(J)
      ENDIF
*
      IPR = J
*
      RETURN
      END
