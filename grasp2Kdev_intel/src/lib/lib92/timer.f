************************************************************************
*                                                                      *
      SUBROUTINE TIMER(I)
*                                                                      *
*   This subroutine  calculates either the cpu time which is left or   *
*   has been used. It has three modes:                                 *
*                                                                      *
*      (1) I = -1  Used to initialise the time stored in  TIME1        *
*      (2) I =  0  Writes out time used in seconds since last call     *
*      (3) I =  1  Returns in  TIME2  the time left in seconds         *
*                                                                      *
*                                          Last updated: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /TIME/TIMMAX,TIMECK,TIME1,TIME2
*
      IF (I .EQ. -1) THEN
*
*   Load  TIME1  with the elapsed CPU time
*
         TIME1 = ZERO
*
      ELSEIF ((I .EQ. 0) .OR. (I .EQ. 1)) THEN
*
*   Load  TIME2  with the elapsed cpu time
*
         TIME2 = ZERO
*
         IF (I .EQ. 0) THEN
*
            TIME2 = TIME2-TIME1
            TIME1 = TIME1+TIME2
*
         ELSEIF (I .EQ. 1) THEN
*
            TIME2 = TEN*TEN
*
         ENDIF
*
      ELSE
*
         WRITE (*,301) I
         STOP
*
      ENDIF
*
      RETURN
*
  301 FORMAT ('TIMER: Argument out of range: I = ',1I8,
     :           '; terminating execution.')
      END
