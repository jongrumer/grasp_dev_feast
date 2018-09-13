************************************************************************
*                                                                      *
      FUNCTION SCREEN (J,UCF)
*                                                                      *
*   This  routine estimates the screening (or shielding) for orbital   *
*   J. The set of rules used by C Froese Fischer in her program MCHF   *
*   are used.                                                          *
*                                                                      *
*   Written by Farid A Parpia, at Oxford  Last revision: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL FULL
*
      DIMENSION UCF(*)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
*
*   Accumulate screening contributions
*
      SCREEN = 0.0D 00
      NPJ = NP(J)
      NLJ = NKL(J)
      DO 1 I = 1,NW
         NPI = NP(I)
         IF (NPJ .EQ. NPI) THEN
            NLI = NKL(I)
            IF (NLJ .GT. NLI) THEN
               SCREEN = SCREEN+UCF(I)
            ENDIF
         ELSEIF (NPJ .GT. NPI) THEN
            SCREEN = SCREEN+UCF(I)
         ENDIF
    1 CONTINUE
*
      UCFJ = UCF(J)
      UCFULL = DBLE (NKJ(J)+1)
      FULL = UCFJ .EQ. UCFULL
      IF (FULL) THEN
         SCREEN = SCREEN+0.5D 00*UCFJ
      ELSEIF (UCFJ .GE. 1.0D 00) THEN
         SCREEN = SCREEN+0.5D 00*(UCFJ-1.0D 00)
      ENDIF
*
      IF (SCREEN .GE. Z) THEN
         WRITE (*,300) SCREEN
         STOP
      ENDIF
*
      RETURN
*
  300 FORMAT ('SCREEN: Screening parameter = ',1P,1D19.12,';'
     :       /' atomic number is probably incorrect')
*
      END
