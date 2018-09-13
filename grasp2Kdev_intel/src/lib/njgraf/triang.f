************************************************************************
*                                                                      *
      SUBROUTINE TRIANG (FAIL)
*                                                                      *
*   Reduces  a triangle having one apex at either end of the axis of   *
*   the flat  diagram.  This introduces one 6j symbol and some phase   *
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
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /NAM/NAMSUB
*
      DATA NAME/'TRIANG'/
*
      NAMSUB = NAME
      IT1 = NPOINT(1)
      IT2 = NPOINT(2)
      IT3 = NPOINT(3)
      JWC = JWC+1
      KW(1,JWC) = JDIAG(IT3,2)
      KW(2,JWC) = JDIAG(IT2,3)
      KW(3,JWC) = JDIAG(IT3,1)
      IF (ARR(IT3,1) .GT. 0) CALL PHASE2 (KW(3,JWC))
      KW(4,JWC) = JDIAG(IT2,1)
      IF (ARR(IT2,1) .LT. 0) CALL PHASE2 (KW(4,JWC))
      K23 = 3
      IF (IT1 .EQ. IFIRST) K23 = 2
      KW(5,JWC) = JDIAG(IT1,K23)
      KW(6,JWC) = JDIAG(IT3,3)
      CALL TRDEL (KW(1,JWC),KW(2,JWC),KW(5,JWC),NBNODE,FAIL)
      IF (FAIL) GOTO 15
      IF (ARR(IT3,3) .GT. 0) CALL PHASE2 (KW(6,JWC))
      JT1 = KW(5,JWC)
      JDIAG(IT3,1) = JT1
      JDIAG(IT3,3) = KW(2,JWC)
      ARR(IT3,1) = ARR(IT1,K23)
      ARR(IT3,3) = ARR(IT2,3)
*
      IF (IT1 .NE. IFIRST) THEN
         TAB1(JT1,1) = IT3
         TAB1(JT1,2) = IH(NBNODE-1)
         K12 = 1
      ELSE
         TAB1(JT1,1) = IH(2)
         TAB1(JT1,2) = IT3
         K12 = 2
      ENDIF
*
      IL3 = IL(IT3)
*
      IF (IT1 .NE. ILAST) THEN
         IL2 = IL(IT2)-1
*
         DO 2 I = 2,IL2
            IT = IH(I)
            ILP = I-1
            IL(IT) = ILP
            IH(ILP) = IT
    2    CONTINUE
      ENDIF
*
      DO 1 I = IL3,NBNODE
         IT = IH(I)
         ILP = I-K12
         IL(IT) = ILP
         IH(ILP) = IT
    1 CONTINUE
*
   15 RETURN
      END
