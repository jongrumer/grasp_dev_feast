************************************************************************
*                                                                      *
      SUBROUTINE CUT1L (FAIL)
*                                                                      *
*   Cut  on  one  line, that  was  left as a free end in JDIAG. Puts   *
*   corresponding delta in J23.                                        *
*                                                                      *
*                                           Last update: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DATA NAME/'CUT1L '/
*
      IT = ITFREE(1)
      J0 = JDIAG(IT,1)
      CALL DELTA (J0,M,FAIL)
      IF (FAIL) GOTO 2
      CALL DELTA (JDIAG(IT,3),JDIAG(IT,2),FAIL)
      IF (FAIL) GOTO 2
      JDIAG(IT+1,3) = JDIAG(IT,3)
*
      IF (ARR(IT,2) .EQ. ARR(IT,3)) THEN
         ARR(IT+1,3) = 1
         ARR(IT-1,2) = -1
      ELSEIF (ARR(IT,2) .LT. ARR(IT,3)) THEN
         ARR(IT+1,3) = -1
         ARR(IT-1,2) = 1
      ENDIF
*
      J9C = J9C+1
      J9(J9C) = JDIAG(IT,3)
      J = 2
      CALL ZERO (J,J0,FAIL)
      IF (FAIL) GOTO 2
      IL1 = IL(IT+1)
*
      DO 1 I = IL1,NBNODE
         IT = IH(I)
         ILP = I-1
         IL(IT) = ILP
         IH(ILP) = IT
    1 CONTINUE
*
      NBNODE = NBNODE-1
*
    2 CALL PRINTJ (NAME,12)
      RETURN
*
      END
