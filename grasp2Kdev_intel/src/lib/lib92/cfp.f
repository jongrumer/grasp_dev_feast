************************************************************************
*                                                                      *
      SUBROUTINE CFP (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,COEFP)
*                                                                      *
*   Selects the appropriate table of  fractional  parentage  coeffi-   *
*   cients in jj-coupling.                                             *
*                                                                      *
*   Input variables:                                                   *
*                                                                      *
*      LOCK     : + OR - (2*j + 1)                                     *
*      NEL      : Number of equivalent electrons in shell              *
*      IJD/IJP  : Total J OF DAUGHTER/PARENT STATE                     *
*      IVD/IVP  : Seniority of daughter/parent state                   *
*      IWD/IWP  : Other quantum number (if needed)                     *
*                                                                      *
*   Output variable:                                                   *
*                                                                      *
*      COEFP    : Numerical result                                     *
*                                                                      *
*   This control routine does not check the input variables for con-   *
*   sistency, except the trivial case  of j = 1/2. All  other checks   *
*   are performed at a lower level. The package will return  correct   *
*   results for j = 3/2, 5/2, 7/2. Higher values of j return a value   *
*   1.0 if NEL = 1 or 2; otherwise 0 with an error signal.             *
*                                                                      *
*   Call(s) to: [LIB92]: CFP3, CFP5, CFP7, CFPD.                       *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
*
      K = IABS (LOCK)/2
*
      IF (K .GE. 5) GOTO 5
*
      GOTO (1,2,3,4),K
*
    1 WRITE (*,300)
      STOP
*
    2 CALL CFP3 (NEL,IJD,IJP,COEFP)
      RETURN
*
    3 CALL CFP5 (NEL,IJD,IVD,IJP,IVP,COEFP)
      RETURN
*
    4 CALL CFP7 (NEL,IJD,IVD,IJP,IVP,COEFP)
      RETURN
*
    5 CALL CFPD (LOCK,NEL,IJD,IVD,IWD,IJP,IVP,IWP,COEFP)
      RETURN
*
  300 FORMAT ('CFP: Unnecessary attempt to form a CFP for an',
     :        ' electron with j = 1/2.')
*
      END
