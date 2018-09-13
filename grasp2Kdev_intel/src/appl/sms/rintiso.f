************************************************************************
*                                                                      *
      FUNCTION RINTISO (I,J)
*                                                                      *
*   The value of RINTISO is an approximation to:                       *
*                                                                      *
*                                                                      *
*         I (dv(r)*( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)        *
*                     I     J       I     J                            *
*                                                                      *
*   where   I ( G(r) ; Range )  denotes  the  integral  of G(r) over   *
*   Range and dv(r) is the difference between the potentials of two    *
*   Fermi nuclear distributions                                        *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*                                                                      *
*   Written by Per Jonsson             Last revision: 24 Dec 1992      *
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
     :      /DVPOT/DV(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Tabulate integrand as required for SUBROUTINE QUAD
*
      MTP = MIN (MF(I),MF(J))
*
*   Value at first tabulation point is arbitrary
*
      TA(1) = 0.0D 00
      DO 1 L = 2,MTP
        TA(L) = DV(L)*(PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))*RP(L)
    1 CONTINUE
*
*   Perform integration
*
      CALL QUAD (RESULT)
      RINTISO = RESULT
*
      RETURN
      END
