************************************************************************
*                                                                      *
      SUBROUTINE CORD (JA,JB,JA1,IPCA,JB1)
*                                                                      *
*   Computes the MCP coefficients for contributions involving closed   *
*   shells.  The  standard formulae are given in I P Grant, Advances   *
*   in Physics  19 (1970) 747, Eq. (8.33).  In this segment JA1, JB1   *
*   point to the JLIST array, IA1, IB1 to the full list of orbitals.   *
*                                                                      *
*   Call(s) to: CLRX, SPEAK.                                           *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
*
      PARAMETER (EPS = 1.0D-10)
*
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
*
*   Set quantum numbers required.
*
      IF (IPCA .EQ. 2) THEN
         IA1 = KLIST(JA1)
      ELSE
         IA1 = JLIST(JA1)
      ENDIF
      IB1 = KLIST(JB1)
*
*   Force IA1 to be greater than IB1
*
      IF (IA1 .GT. IB1) THEN
         NS = IA1
         IA1 = IB1
         IB1 = NS
      ENDIF
*
      KAP1 = NAK(IA1)
      J1 = IABS (KAP1)
      NQS1 = NQ1(IA1)
*
      IF (IA1 .EQ. IB1) THEN
*
*   Case when IA1 .EQ. IB1
*
         X = DBLE (NQS1*(NQS1-1)/2)
         CALL SPEAK (JA,JB,IA1,IB1,IA1,IB1,0,X)
         NUMAX = J1+J1-2
         IF (NUMAX  .LE.  0) RETURN
         CONST = DBLE (NQS1*NQS1/2)
         DO 1 NU = 2,NUMAX,2
            GAM = CLRX (KAP1,NU,KAP1)
            X = -CONST*GAM*GAM
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IA1,IB1,NU,X)
    1    CONTINUE
*
*   Case when IA1 .NE. IB1
*
      ELSE
*
         KAP2 = NAK(IB1)
         J2 = ABS (KAP2)
         NQS2 = NQ1(IB1)
         CONST = DBLE (NQS1*NQS2)
         CALL SPEAK (JA,JB,IA1,IB1,IA1,IB1,0,CONST)
         NUMIN = ABS (J1-J2)
         NUMAX = J1+J2-1
         IF (KAP1*KAP2 .LT. 0) NUMIN = NUMIN+1
         DO 2 NU = NUMIN,NUMAX,2
            GAM = CLRX (KAP1,NU,KAP2)
            X = -CONST*GAM*GAM
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IB1,IA1,NU,X)
    2    CONTINUE
*
      ENDIF
*
      RETURN
      END
