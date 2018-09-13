************************************************************************
*                                                                      *
      SUBROUTINE SETMC
*                                                                      *
*   This subprogram sets machine-dependent parameters.                 *
*                                                                      *
*   Call(s) to: [LAPACK]: DLAMCH.                                      *
*                                                                      *
*   Written by Farid A Parpia              Last updated: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL LDBPG
*
      COMMON/DEBUGG/LDBPG(5)
     :      /DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
*
*   Set the machine-dependent parameters:
*
*   TENMAX - Maximum size of exponent of 10
*
      TENMAX = DLAMCH ('L')
*
*   EXPMAX - Maximum size of exponent of e
*
      DNUM = DLAMCH ('O')
      EXPMAX = LOG (DNUM)
*
*   EXPMIN - Minimum size of exponent of e
*
      DNUM = DLAMCH ('U')
      EXPMIN = LOG (DNUM)
*
*   PRECIS - Machine precision
*
      PRECIS = DLAMCH ('E')
*
*   Debug printout
*
      IF (LDBPG(1)) WRITE (99,300) TENMAX,EXPMAX,EXPMIN,PRECIS
*
      RETURN
*
  300 FORMAT (/'From SUBROUTINE SETMC:'
     :        /' TENMAX (maximum exponent of 10): ',F5.0
     :        /' EXPMAX (maximum exponent of e): ',1P,1D19.12
     :        /' EXPMIN (minimum exponent of e): ',   1D19.12
     :        /' PRECIS (machine precision): ',1D19.12)
*
      END
