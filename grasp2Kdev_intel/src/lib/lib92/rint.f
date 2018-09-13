************************************************************************
*                                                                      *
      FUNCTION RINT (I,J,K)
*                                                                      *
*   The value of RINT is an approximation to:                          *
*                                                                      *
*              k                                                       *
*         I ( r  *  ( P (r)*P (r) + Q (r)*Q (r) ; 0 to infinity)       *
*                      i     j       i     j                           *
*                                                                      *
*   where   I ( G(r) ; range )  denotes  the  integral  of G(r) over   *
*   range.                                                             *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 05 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Determine the maximum tabulation point for the integrand
*
      MTP = MIN (MF(I),MF(J))
*
*   Tabulate the integrand as required for SUBROUTINE QUAD; the
*   value at the first tabulation point is arbitrary
*
      TA(1) = 0.0D 00
      DO 1 L = 2,MTP
         TA(L) = (R(L)**K)*(PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))*RP(L)
    1 CONTINUE
*
*   Perform the quadrature
*
      CALL QUAD (RESULT)
      RINT = RESULT
*
      RETURN
*
      END
