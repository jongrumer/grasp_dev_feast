************************************************************************
*                                                                      *
      FUNCTION RINT_SMS2 (I,J)
*                                                                      *
*   The value of RINT_SMS2 is an approximation to:                   *
*                                                                      *
*                                                                      *
*   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
*   Range.                                                             *
*                                                                      *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*                                                                      *
*   Written by  G. Gaigalas                                            *
*           and E. Gaudamauskas                                        *
*                        Last             revision:  09 October 2009   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
*
      POINTER (PNTRPF,PF(NNNP,*)),(PNTRQF,QF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /NPAR/PARM(2),NPARM
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF2/C
CGG
      COMMON/DEBUGG/LDBPG(5)
     :      /DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
*
*   Tabulate integrand as required for SUBROUTINE QUAD
*
      MTP = MIN (MF(I),MF(J))
*
*   Value at first tabulation point is arbitrary
*
      CALL SIGMA_1(2,I,J,APART1)
      CALL SIGMA_2(2,I,J,APART2)
      TA(1) = 0.0D 00
      DO 1 L = 2,MTP
*HFS        TA(L) = (R(L)**K)*(PF(L,I)*QF(L,J)+QF(L,I)*PF(L,J))*RP(L)
        TA(L) =  RP(L)*(APART1*QF(L,I)*PF(L,J)
     :        -         APART2*PF(L,I)*QF(L,J))/R(L)
    1 CONTINUE
*
*   Perform integration
*
      CALL QUAD (RESULT)
      RINT_SMS2 = -RESULT*Z/C
*
      RETURN
      END
