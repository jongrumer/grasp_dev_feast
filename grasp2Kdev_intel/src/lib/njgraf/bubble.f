************************************************************************
*                                                                      *
      SUBROUTINE BUBBLE (JPOL,FAIL)
*                                                                      *
*   Reduces a circuit of order 2 , giving  delta  function and phase   *
*   factors.                                                           *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME,NAMSUB
      INTEGER ARR,TAB1
      LOGICAL FAIL,SUMVAR
*
      PARAMETER (
     :   MANGM = 60, M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/NAM/NAMSUB
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'BUBBLE'/
*
      NAMSUB = NAME
      K2 = 2
      K23 = 3
      I1 = 1
      I2 = 1
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
*
      IF (IT2 .EQ. ILAST) THEN
         IF (IT1 .NE. IFIRST) THEN
            IT2 = IT1
            IT1 = ILAST
         ENDIF
         I1 = -1
         K23 = 2
         I2 = 2
      ENDIF
*
      CALL PHASE (IT1,JDIAG,M4TRD)
      K = ABS ((3*ARR(IT2,1)+2*ARR(IT2,2)+ARR(IT2,3))/2)+1
      IF (K .NE. 4) CALL PHASE2 (JDIAG(IT2,K))
      IF (NBNODE .EQ. 2) RETURN
      IL1 = IL(IT2)+I1
      IT = IH(IL1)
      ARR(IT,K23) = ARR(IT1,K23)
      L = JDIAG(IT1,K23)
      L1 = JDIAG(IT,K23)
      JDIAG(IT,K23) = L
*
      IF (JPOL .NE. 1) THEN
         CALL DELTA (L,L1,FAIL)
         IF (FAIL) RETURN
      ELSE
         MP = MP-1
         KW(2,JWC) = L
         J6(J6C-1) = L
         J6(J6C) = L
         IF (K .EQ. 2) J8(J8C) = L
      ENDIF
*
      TAB1(L,I2) = IT
*
      IF (IT1 .NE. ILAST) THEN
         IF (IT2 .EQ. ILAST) THEN
            TAB1(L,1) = IH(2)
            IL1 = 2
            K2 = 1
         ENDIF
*
         DO 1 I = IL1,NBNODE
            IT = IH(I)
            IL(IT) = I-K2
            IH(I-K2) = IT
    1    CONTINUE
*
      ENDIF
*
      J9(J9C+1) = L
      J9C = J9C+2
      J9(J9C) = L
*
      RETURN
      END
