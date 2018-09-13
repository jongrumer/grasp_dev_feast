************************************************************************
*                                                                      *
      FUNCTION ARCTAN (ARG1,ARG2)
*                                                                      *
*                   -1                                                 *
*       ARCTAN = tan   (ARG1/ARG2),      0 .LE. ARCTAN .LT. 2*\pi      *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL FIRST,INTRIN
*
      COMMON/DEF9/CVAC,PI
*
      DATA FIRST/.TRUE./,
     :     INTRIN/.TRUE./
*
*   Determine whether the FORTRAN intrinsic function ATAN2 always
*   returns a positive value
*
      IF (FIRST) THEN
         ARCTAN = ATAN2 (-1.0D 00,-1.0D 00)
         IF (ARCTAN .GT. 0.0D 00) THEN
            INTRIN = .TRUE.
         ELSE
            INTRIN = .FALSE.
         ENDIF
         FIRST = .FALSE.
      ENDIF
*
*   Use the intrinsic function if it passes the above test; otherwise
*   add 2*PI to the negative values returned by the intrinsic function
*
      IF (INTRIN) THEN
         ARCTAN = ATAN2 (ARG1,ARG2)
      ELSE
         IF (ARG1 .GE. 0.0D 00) THEN
            ARCTAN = ATAN2 (ARG1,ARG2)
         ELSE
            ARCTAN = PI+PI+ATAN2 (ARG1,ARG2)
         ENDIF
      ENDIF
*
      RETURN
      END
