************************************************************************
*                                                                      *
      FUNCTION SLATER (IA,IB,IC,ID,K)
*                                                                      *
*   The value of this  function is the Slater integral                 *
*                                                                      *
*                               k                                      *
*                              R (abcd)                                *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD, YZK.                                    *
*                                                                      *
*                                         Last revision: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      CALL YZK (K,IB,ID)
*
*   Multiply by second term, and obtain result by integration
*
      IF ((K .EQ. 0) .AND. (IB .EQ. ID)) THEN
         MTP = MIN (    MF(IA),MF(IC))
      ELSE
         MTP = MIN (MTP,MF(IA),MF(IC))
      ENDIF
*
      TA(1) = 0.0D 00
      DO 1 I = 2,MTP
         TA(I) =  (PF(I,IA)*PF(I,IC)+QF(I,IA)*QF(I,IC))
     :           *RPOR(I)*TB(I)
    1 CONTINUE
*
      CALL QUAD (RESULT)
      SLATER = RESULT
*
*   Debug printout
*
      IF (LDBPR(10)) WRITE (99,300) K,NP(IA),NH(IA),NP(IB),NH(IB),
     :                               NP(IC),NH(IC),NP(ID),NH(ID),
     :                                 SLATER
*
      RETURN
*
  300 FORMAT (/'  (',1I1,')',
     :        /' R   (',1I2,1A2,',',1I2,1A2,';',1I2,1A2,',',1I2,1A2,') '
     :        ,'= ',1PD19.12)
*
      END
