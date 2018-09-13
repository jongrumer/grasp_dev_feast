************************************************************************
*                                                                      *
      SUBROUTINE INTAB
*                                                                      *
*   This SUBROUTINE called at the end of DIAGRM, fixes the arrays IH   *
*   and IL - so to speak hardware and logical addresses of triads in   *
*   JDIAG . Also  determines the number of free ends NFREE and their   *
*   location ITFREE.                                                   *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INTEGER ARR,TAB1
      LOGICAL FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DO 1 I = 1,M
         IAL(I) = 1
    1 CONTINUE
*
      DO 3 I = IFIRST,ILAST
         J = JDIAG(I,1)
         K = IAL(J)
         TAB1(J,K) = I
         IAL(J) = K+1
    3 CONTINUE
*
      IFR = IFIRST-1
*
      DO 4 I = IFIRST,ILAST
         IT = I-IFR
         IL(I) = IT
         IH(IT) = I
    4 CONTINUE
*
      J = JDIAG(IFIRST,3)
      K = IAL(J)
      IF (K .GT. 1) TAB1(J,2) = TAB1(J,1)
      TAB1(J,1) = IFIRST
      IAL(J) = 3
      J = JDIAG(ILAST,2)
      TAB1(J,2) = ILAST
      IAL(J) = 3
      NFREE = 0
*
      DO 7 I = IFIRST,ILAST
         J = JDIAG(I,1)
         IF (IAL(J) .NE. 3) THEN
            NFREE = NFREE+1
            ITT = ILAST+NFREE
            TAB1(J,2) = ITT
            IL(ITT) = NFREE*1000
            ITFREE(NFREE) = I
         ENDIF
    7 CONTINUE
*
      RETURN
      END
