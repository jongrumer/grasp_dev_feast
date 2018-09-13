************************************************************************
*                                                                      *
      SUBROUTINE LOLPOP (FAIL)
*                                                                      *
*   Reduces a loop with one line and one node in the flat graph.       *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION KP(3),KS(3)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /NAM/NAMSUB
*
      DATA NAME/'LOLPOP'/
      DATA KP/2,3,1/
      DATA KS/0,1,-1/
*
      NAMSUB = NAME
      I1 = NPOINT(1)
      K3 = 2
      IF (I1 .EQ. ILAST) K3 = 3
      L = JDIAG(I1,K3)
      CALL DELTA (L,MP,FAIL)
      IF (FAIL) RETURN
      K = KP(K3)
      IF (ARR(I1,K) .LT. 0) CALL PHASE2 (JDIAG(I1,K))
      K1 = KS(K3)
      IL1 = IL(I1)+K1
      I2 = IH(IL1)
      L1 = JDIAG(I2,1)
      CALL DELTA (L1,JDIAG(I2,K3),FAIL)
      IF (FAIL) RETURN
      IF (ARR(I2,K3) .EQ. K1) CALL PHASE2 (L1)
      IL2 = IL(I2)+K1
      I3 = IH(IL2)
      K2 = K3+K1
      JDIAG(I3,K2) = L1
      ARR(I3,K2) = ARR(I2,1)
      J9C = J9C+1
      J9(J9C) = L1
      J6C = J6C+1
      J6(J6C) = JDIAG(I1,1)
      IF (K3 .EQ. 3) RETURN
*
      DO 1 I = 3,NBNODE
         IT = IH(I)
         ILP = I-2
         IL(IT) = ILP
         IH(ILP) = IT
    1 CONTINUE
*
      RETURN
      END
