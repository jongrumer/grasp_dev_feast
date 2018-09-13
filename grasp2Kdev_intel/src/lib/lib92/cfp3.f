************************************************************************
*                                                                      *
      SUBROUTINE CFP3 (NEL,IJD,IJP,COEFP)
*                                                                      *
*   Table look-up for fractional parentage coefficients of  equival-   *
*   ent electrons with j = 3/2. See listing of CFP for argument list.  *
*                                                                      *
*                                           Last update: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
*   Floating point constants
*
      PARAMETER (C1 = 6.0D 00,C2 = 5.0D 00)
*
      IF ((NEL .LE. 0) .OR. (NEL .GT. 4)) GOTO 7
*
      GOTO (1,2,3,5), NEL
    1 IF ((IJD .NE. 3) .OR.  (IJP .NE. 0)) GOTO 7
      GOTO 6
*
    2 IF (IJP .NE. 3) GOTO 7
      IF ((IJD .EQ. 0) .OR. (IJD .EQ. 4)) GOTO 6
      GOTO 7
*
    3 IF (IJD .NE. 3) GOTO 7
      IF (IJP .NE. 0) GOTO 4
      COEFP = SQRT (1.0D 00/C1)
      RETURN
    4 IF (IJP .NE. 4) GOTO 7
      COEFP = -SQRT (C2/C1)
      RETURN
*
    5 IF ((IJD .NE. 0) .OR. (IJP .NE. 3)) GOTO 7
    6 COEFP = 1.0D 00
      RETURN
*
*   Fault mode section
*
    7 WRITE (*,300) NEL,IJD,IJP
      STOP
*
  300 FORMAT ('CFP3: Error in trying to compute a CFP',
     :        ' for a state with ',1I2,' electrons with j = 3/2;'
     :       /' IJD = ',1I2,', IJP = ',1I2,'.')
*
      END
