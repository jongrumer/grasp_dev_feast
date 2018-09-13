************************************************************************
*                                                                      *
      FUNCTION RINTII (I,J,K)
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
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)

      include 'parameters.def'
CGG      PARAMETER (NNNP = 590) 
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      POINTER (PNTRPFII,PFII(NNNP,*)),(PNTRQFII,QFII(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVEII/PZII(NNNW),PNTRPFII,PNTRQFII,MFII(NNNW)
*
*   Determine the maximum tabulation point for the integrand
*
      MTP = MIN (MFII(I),MFII(J))
*
*   Tabulate the integrand as required for SUBROUTINE QUAD; the
*   value at the first tabulation point is arbitrary
*
      TA(1) = 0.0D 00
      DO L = 2,MTP
         TA(L) = (R(L)**K)*(PFII(L,I)*PFII(L,J)+
     :           QFII(L,I)*QFII(L,J))*RP(L)
      ENDDO
*
*   Perform the quadrature
*
      CALL QUAD (RESULT)
      RINTII = RESULT
*
      RETURN
*
      END
