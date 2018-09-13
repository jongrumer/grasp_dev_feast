************************************************************************
*                                                                      *
      FUNCTION OCON (IA1,IB1,IA2,IB2)
*                                                                      *
*   Evaluates the  multiplicative statistical  factor. It is assumed   *
*   that states are ordered so that IA1 .LE. IB1, IA2 .LE. IB2.        *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
*
      COMMON/M1/NQ1(NNNW),NQ2(NNNW)
*
      WA = DBLE (NQ1(IA1)*NQ1(IB1))
      IF (IA1 .EQ. IB1) WA = WA-DBLE (NQ1(IA1))
      WB = DBLE (NQ2(IA2)*NQ2(IB2))
      IF (IA2 .EQ. IB2) WB = WB-DBLE (NQ2(IB2))
      WC = WA*WB
      OCON = SQRT (WC)
*
*   Set phase factor (-1)**(DELTA P)
*
      LRD1 = MIN (IA2,IB2)+1
      LRD2 = MAX (IA2,IB2)
      IF (LRD1 .GT. LRD2) THEN
         IDR = 0
      ELSE
         IDR = 1
         DO 1 K = LRD1,LRD2
            IDR = IDR+NQ2(K)
    1    CONTINUE
      ENDIF
*
      LLD1 = MIN (IA1,IB1)+1
      LLD2 = MAX (IA1,IB1)
      IF (LLD1 .GT. LLD2) THEN
         IDL = 0
      ELSE
         IDL = 1
         DO 2 K = LLD1,LLD2
            IDL = IDL+NQ1(K)
    2    CONTINUE
      ENDIF
*
      IPHAS = IDR-IDL
      IF (MOD (IPHAS,2) .NE. 0) OCON = -OCON
*
      RETURN
      END
