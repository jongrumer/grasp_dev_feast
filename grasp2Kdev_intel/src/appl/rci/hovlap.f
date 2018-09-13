************************************************************************
*                                                                      *
      FUNCTION HOVLAP (P,Q,MTPO,NP,KAPPA,Z)
*                                                                      *
*   This subprogram computes the overlap of the orbital tabulated in   *
*   the arrays  P  and  Q  with maximum tabulation point  MTPO  with   *
*   a hydrogenic orbital with parameters  NP  KAPPA  Z .               *
*                                                                      *
*   Call(s) to: [LIB92]: DCBSRW, QUAD.                                 *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      DIMENSION P(NNNP),Q(NNNP)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /HORB/PH(NNNP),QH(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Set up the hydrogenic orbital
*
      CALL DCBSRW (NP,KAPPA,Z,EH,PZH,PH,QH,MTPH)
*
*   Compute the overlap
*
      MTP = MIN (MTPH,MTPO)
      TA(1) = 0.0D 00
      DO 1 I = 2,MTP
          TA(I) = (P(I)*PH(I)+Q(I)*QH(I))*RP(I)
    1 CONTINUE
      CALL QUAD (RESULT)
*
      HOVLAP = RESULT
*
      RETURN
      END
