************************************************************************
*                                                                      *
      FUNCTION RINTFF (I,J,K)
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
      POINTER (PNTRPFFF,PFFF(NNNP,*)),(PNTRQFFF,QFFF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVEFF/PZFF(NNNW),PNTRPFFF,PNTRQFFF,MFFF(NNNW)
*
*   Determine the maximum tabulation point for the integrand
*
      MTP = MIN (MFFF(I),MFFF(J))
*
*   Tabulate the integrand as required for SUBROUTINE QUAD; the
*   value at the first tabulation point is arbitrary
*
      TA(1) = 0.0D 00
      DO L = 2,MTP
         TA(L) = (R(L)**K)*(PFFF(L,I)*PFFF(L,J)+
     :           QFFF(L,I)*QFFF(L,J))*RP(L)
      ENDDO
*
*   Perform the quadrature
*
      CALL QUAD (RESULT)
      RINTFF = RESULT
*
      RETURN
*
      END
