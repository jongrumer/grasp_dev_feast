************************************************************************
*                                                                      *
      SUBROUTINE CUTNL (FAIL)
*                                                                      *
*   This subroutine  examines the case where there are more than two   *
*   free ends, but they are contiguous, so that the graph can be cut   *
*   without destroying the flat structure.                             *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW,ARR,TAB1
      LOGICAL TABS,SUMVAR,FAIL
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'CUTNL '/
*
      NTF = ITFREE(NFREE)-ITFREE(1)
      IF (NTF .GT. NFREE) GOTO 8
      IT2 = ITFREE(1)
      IT3 = ITFREE(NFREE)
      IT1 = IT2-1
      IT4 = IT3+1
*
      IF (NTF .NE. NFREE) THEN
*
         JT = JDIAG(IT2,3)
         CALL DELTA (JT,JDIAG(IT3,2),FAIL)
*
         IF (FAIL) GOTO 9
*
         IF (ARR(IT2,3) .EQ. ARR(IT3,2)) THEN
            CALL PHASE2 (JT)
            ARR(IT2,3) = -ARR(IT2,3)
            ARR(IT1,2) = -ARR(IT1,2)
         ENDIF
*
         JDIAG(IT3,2) = JT
         JDIAG(IT4,3) = JT
         J9(J9C+1) = JT
         J9C = J9C+2
         J9(J9C) = JT
         NBTR = NBTR+NFREE
         IT5 = 0
*
      ELSE
*
         NFR = 0
*
         DO 3 IT5 = IT2,IT3
            NFR = NFR+1
            IF (ITFREE(NFR) .GT. IT5) GOTO 4
    3    CONTINUE
*
    4    JKP(1,1) = JDIAG(IT5,1)
         JARR(1,1) = -ARR(IT5,1)
         JKP(1,2) = JDIAG(IT2,3)
         JARR(1,2) = -ARR(IT2,3)
         JKP(1,3) = JDIAG(IT3,2)
         JARR(1,3) = -ARR(IT3,2)
*
         DO 5 J = 1,3
            JKP(2,J) = JDIAG(IT5,J)
            JARR(2,J) = ARR(IT5,J)
    5    CONTINUE
*
         JDIAG(IT5,2) = JDIAG(IT3,2)
         ARR(IT5,2) = ARR(IT3,2)
         JDIAG(IT5,3) = JDIAG(IT2,3)
         ARR(IT5,3) = ARR(IT2,3)
         ILP = IL(IT2)
         IL(IT5) = ILP
         IH(ILP) = IT5
         NBTR = NBTR+NFREE+2
         CALL PHASE (IT5,JDIAG,M4TRD)
Cww               
Cww         ERROR in NJGRAF
Cww         K = ABS ((3*ARR(IT5,1)+2*ARR(IT5,2)+ARR(IT5,3))/2+1)
Cww               
         K = ABS (3*ARR(IT5,1)+2*ARR(IT5,2)+ARR(IT5,3))/2+1
         IF (K .NE. 4) CALL PHASE2 (JDIAG(IT5,K))
*
      ENDIF
*
      IL1 = IL(IT4)
*
      DO 7 I = IL1,NBNODE
         IT = IH(I)
         ILP = I-NFREE
         IL(IT) = ILP
         IH(ILP) = IT
    7 CONTINUE
*
      NBNODE = NBNODE-NFREE
      NFIN = 0
      GOTO 8
*
    9 FAIL = .TRUE.
    8 CALL PRINTJ (NAME,8)
*
      RETURN
*
      END
