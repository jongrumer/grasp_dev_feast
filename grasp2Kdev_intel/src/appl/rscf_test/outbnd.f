************************************************************************
*                                                                      *
      FUNCTION OUTBND (ETRY)
*                                                                      *
*   This  subprogram determines whether the trial eigenvalue etry is   *
*   within the bounds (EPSMIN,EPSMAX)                                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL OUTBND
*
      COMMON/SCF4/EPSMIN,EPSMAX,EMIN,EMAX,ZINF
*
      IF ((ETRY .GT. EPSMIN) .AND. (ETRY .LT. EPSMAX)) THEN
         OUTBND = .FALSE.
      ELSE
         OUTBND = .TRUE.
      ENDIF
*
      RETURN
      END
