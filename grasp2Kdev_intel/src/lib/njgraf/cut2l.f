************************************************************************
*                                                                      *
      SUBROUTINE CUT2L (FAIL)
*                                                                      *
*   Cut on two lines that were left as free ends in JDIAG. Puts cor-   *
*   responding delta in J23.                                           *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1,ARROW
      LOGICAL FAIL,TABS,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,    MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'CUT2L '/
*
*
      IT1 = ITFREE(1)
      IT2 = ITFREE(2)
      JT1 = JDIAG(IT1,1)
      JT2 = JDIAG(IT2,1)
      CALL DELTA (JT1,JT2,FAIL)
      IF (FAIL) GOTO 1
      IF (ARR(IT1,1) .EQ. ARR(IT2,1)) CALL PHASE2 (JT1)
      ARR(IT2,1) = -ARR(IT1,1)
      JDIAG(IT2,1) = JT1
      TAB1(JT1,2) = IT2
      J9(J9C+1) = JT1
      J9C = J9C+2
      J9(J9C) = JT1
      CALL OTHERJ (0,JT1,L1,LC1,K1)
      CALL OTHERJ (0,JT2,L2,LC2,K2)
      J23(L2,LC2) = JT1
      LINE(JT1,K1) = L2
      LCOL(JT1,K1) = LC2
      ARROW(L2,LC2) = -ARROW(L1,LC1)
*
    1 CALL PRINTJ (NAME,12)
*
      RETURN
      END
