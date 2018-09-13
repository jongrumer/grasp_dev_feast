************************************************************************
*                                                                      *
      SUBROUTINE CFPD (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,COEFP)
*                                                                      *
*   This is a dummy SUBROUTINE. It returns correct values for 1 or 2   *
*   particle or single hole states, and signals an error otherwise.    *
*                                                                      *
*                                           Last update: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
      IF ((NEL .EQ. 1)
     :    .OR. (NEL .EQ. 2)
     :      .OR. (NEL .EQ. ABS (LOCK))) THEN
         COEFP = 1.0D 00
      ELSE
         LOCJ = IABS (LOCK) - 1
         WRITE (*,300) LOCJ
         STOP
      ENDIF
*
      RETURN
*
  300 FORMAT ('CFPD: Inadmissable attempt to obtain a CFP',
     :        ' for a state of a shell with j = ',1I2,'/2;'
     :       /' new subprogram required.')
*
      END
