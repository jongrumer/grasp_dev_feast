************************************************************************
*                                                                      *
      SUBROUTINE DEFCOR (J)
*                                                                      *
*   Compute the deferred corrections for orbital J .                   *
*                                                                      *
*                                          Last updated: 18 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL FIRST
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /DEFC/DP(NNNP),DQ(NNNP)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      PARAMETER (W3 =    1.0D 00/120.0D 00,
     :           W2 =  -15.0D 00*W3,
     :           W1 =   40.0D 00*W3)
*
      DATA FIRST /.TRUE./
*
*   The deferred corrections for the first two points are
*   unnecessary, because the integration always commences
*   from the fourth point of the grid; this requires the
*   deferred correction at the third and subsequent points
*   only
*
      IF (FIRST) THEN
         DO 1 I = 1,2
            DP(I) = 0.0D 00
            DQ(I) = 0.0D 00
    1    CONTINUE
         FIRST = .FALSE.
      ENDIF
*
*   Intermediate points
*
      MFJM3 = MF(J)-3
      DO 2 I = 3,MFJM3
*
         DP(I) = W3*(PF(I+3,J)-PF(I-2,J))
     :          +W2*(PF(I+2,J)-PF(I-1,J))
     :          +W1*(PF(I+1,J)-PF(I  ,J))
*
         DQ(I) = W3*(QF(I+3,J)-QF(I-2,J))
     :          +W2*(QF(I+2,J)-QF(I-1,J))
     :          +W1*(QF(I+1,J)-QF(I  ,J))
*
    2 CONTINUE
*
*   Set remaining deferred corrections to zero: slopes are
*   small in this region
*
      MFJM2 = MF(J)-2
      DO 3 I = MFJM2,N
         DP(I) = 0.0D 00
         DQ(I) = 0.0D 00
    3 CONTINUE
*
      RETURN
      END
