************************************************************************
*                                                                      *
      SUBROUTINE GETCPR (RRMS,APARM,CPARM)
*                                                                      *
*   Determines the parameter `c' (CPARM) for a Fermi nucleus,  given   *
*   the root mean square radius (RRMS) and the parameter `a' (APARM).  *
*   We use the formalism developed in F. A. Parpia and A. K. Mohanty   *
*   ``Relativistic basis set  calculations for atoms with  Fermi nu-   *
*   clei'' Phys Rev A (1992) in press.                                 *
*                                                                      *
*   Call(s) to: ESTRMS.                                                *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
cxhb
cxh Accuracy parameter
cxhe
      ACCY = 1.0D-12
*
*   Bracket CPARM with a lower and upper limit
*
*   Lower limit
*
      CPMIN = RRMS
    1 CPMIN = 0.5D 00*CPMIN
      IF (ESTRMS (APARM,CPMIN) .GT. RRMS) GOTO 1
*
*   Upper limit
*
      CPMAX = RRMS
    2 CPMAX = 2.0D 00*CPMAX
      IF (ESTRMS (APARM,CPMAX) .LT. RRMS) GOTO 2
*
*   Find CPARM by the method of bisection
*
    3 CPTRY = 0.5D 00*(CPMAX+CPMIN)
*
      RMSTRY = ESTRMS (APARM,CPTRY)
*
      IF (RMSTRY .GT. RRMS) THEN
         CPMAX = CPTRY
      ELSE
         CPMIN = CPTRY
      ENDIF
*
      IF ( ( (CPMAX-CPMIN)/(CPMAX+CPMIN) .GT. ACCY ) .AND.
     :     ( ABS (RMSTRY-RRMS)/RRMS      .GT. ACCY ) ) GOTO 3
*
      CPARM = CPTRY
*
      RETURN
      END
