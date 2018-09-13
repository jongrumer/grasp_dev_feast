************************************************************************
*                                                                      *
      SUBROUTINE DPBDT (j)
*                                                                      *
*   This subroutine computes H times the derivative, with respect to   *
*   the internal grid, of the large and small components of the wave   *
*   function with index  J .  These  are tabulated, respectively, in   *
*   arrays  TA  and  TB  in  COMMON  block  /TATB/ .                   *
*                                                                      *
*   A  thirteen-point  Lagrange formaula is used for the calculation   *
*   of derivatives.                                                    *
*                                                                      *
*   Written by Farid F Parpia, at Oxford   Last updated: 06 Oct 1992   *
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
     :      /LIC13/A(13,13)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      EQUIVALENCE (A1,A(7,1)),(A2,A(7,2)),(A3,A(7,3)),
     :            (A4,A(7,4)),(A5,A(7,5)),(A6,A(7,6))
*
*   Compute derivative in three separate regions
*
*   First, points 1 to 6
*
      DO 2 I = 1,6
         HDPBDT = 0.D0
         HDQBDT = 0.D0
         DO 1 K = 1,13
            AIK = A(I,K)
            HDPBDT = HDPBDT+AIK*PF(K,J)
            HDQBDT = HDQBDT+AIK*QF(K,J)
    1    CONTINUE
         TA(I) = HDPBDT
         TB(I) = HDQBDT
    2 CONTINUE
*
*   Next, points 7 to N-6
*
*   Special treatment for this region because of the symmetry of
*   the differentiation formula
*
      DO 3 I = 7,N-6
         TA(I) =  A1*(PF(I-6,J)-PF(I+6,J))+A2*(PF(I-5,J)-PF(I+5,J))
     :           +A3*(PF(I-4,J)-PF(I+4,J))+A4*(PF(I-3,J)-PF(I+3,J))
     :           +A5*(PF(I-2,J)-PF(I+2,J))+A6*(PF(I-1,J)-PF(I+1,J))
         TB(I) =  A1*(QF(I-6,J)-QF(I+6,J))+A2*(QF(I-5,J)-QF(I+5,J))
     :           +A3*(QF(I-4,J)-QF(I+4,J))+A4*(QF(I-3,J)-QF(I+3,J))
     :           +A5*(QF(I-2,J)-QF(I+2,J))+A6*(QF(I-1,J)-QF(I+1,J))
    3 CONTINUE
*
*   Last, points N-5 to N
*
      DO 5 I = N-5,N
         IROW = I-N+13
         HDPBDT = 0.D0
         HDQBDT = 0.D0
         DO 4 K = 1,13
            AIK = A(IROW,K)
            LOC = N-13+K
            HDPBDT = HDPBDT+AIK*PF(LOC,J)
            HDQBDT = HDQBDT+AIK*QF(LOC,J)
    4    CONTINUE
         TA(I) = HDPBDT
         TB(I) = HDQBDT
    5 CONTINUE

      RETURN

      END
