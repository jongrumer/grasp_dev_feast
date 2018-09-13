************************************************************************
*                                                                      *
      SUBROUTINE CHANGE (L,K)
*                                                                      *
*   Exchanges the free ends in either first or last triad of JDIAG.    *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INTEGER ARR,TAB1
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM,
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,    MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
*
      CALL PHASE (L,JDIAG,M4TRD)
      JP = JDIAG(L,K)
      JDIAG(L,K) = JDIAG(L,1)
      JDIAG(L,1) = JP
      JAR = ARR(L,K)
      ARR(L,K) = ARR(L,1)
      ARR(L,1) = JAR
*
      RETURN
      END
