************************************************************************
*                                                                      *
      SUBROUTINE ORDTRI
*                                                                      *
*   This subroutine orders the triads which were left with free ends   *
*   as consequence of cutting,so that the new graph will start there.  *
*                                                                      *
*                                           Last update: 16 Aug 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW,ARR,TAB1
      LOGICAL SUMVAR,TABS
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /GRAPH/JDIAG(M4TRD,3),ARR(M4TRD,3),TAB1(MANGM,2),IL(M4TRD),
     :             IH(M4TRD),NPOINT(M2TRD),NBNODE,IFIRST,ILAST,IPARTS,
     :             IPARTL,NPART,ICROSS,NFREE,ITFREE(M6J),NFIN,NC
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /KEEP/JKP(2,3),JARR(2,3),IT2,IT3,IT5
*
      DATA NAME/'ORDTRI'/
*
      DO 10 I = 1,MP
         IAL(I) = 0
   10 CONTINUE
*
      IF (NFIN .NE. 0) THEN
         NBT1 = NBTR-1
         NBT = NBT1+NFIN
         NBTT = NBT+1
         NB = 0
         GOTO 31
      ENDIF
*
      NF = NBTR-ITFREE(1)
*
      IF (IT5 .EQ. 0) THEN
         NBT1 = NBTR-1
         N0 = 0
         NFT = NFREE
         ISW = 2
         GOTO 100
      ENDIF
*
      NFT = IT5-IT2
      NM = NFT+NBTR+1
      NBT1 = NBTR
*
      DO 21 J = 1,3
         JDIAG(NBTR,J) = JKP(1,J)
         ARR(NBTR,J) = JARR(1,J)
   21 CONTINUE
*
      JT = JDIAG(NM,1)
      N0 = 0
      ISW = 1
      GOTO 100
*
   22 N0 = NFT
*
      DO 211 J = 1,3
         JDIAG(NM,J) = JKP(2,J)
         ARR(NM,J) = JARR(2,J)
  211 CONTINUE
*
      NBT1 = NBT1+1
      NFT = IT3-IT5
      ISW = 3
      GOTO 100
*
   24 NBT1 = K-NFT
*
   23 NODE = NBT1+NFT
      CALL CHANGE (NODE,2)
      GOTO 40
*
   31 DO 35 I = 1,NBNODE
         I1 = IH(I)
         IF (IL(I1) .GT. ILAST) GOTO 35
         I2 = NBT1+I
         IF (I1 .GT. NBTT) GOTO 33
         IF (I1 .EQ. I2) GOTO 32
         IF (IL(I2) .LE. NBNODE) GOTO 35
*
   33    DO 34 J = 1,3
            JDIAG(I2,J) = JDIAG(I1,J)
            ARR(I2,J) = ARR(I1,J)
   34    CONTINUE
*
         IL(I1) = ILAST+I
   32    NB = NB+1
         IL(I2) = 0
*
   35 CONTINUE
*
      IF (NB .NE. NFIN) GOTO 31
      NODE = NBT
   40 IF1 = JDIAG(NBTR,1)
      IF2 = JDIAG(NBTR,3)
*
      DO 51 I = NBTR,NODE
         DO 50 K = 1,3
            J = JDIAG(I,K)
            IAL(J) = IAL(J)+1
   50    CONTINUE
   51 CONTINUE
*
      ILAST = NODE
      CALL PRINTJ (NAME,8)
*
      RETURN
*
  100 IF (NF .LE. 0) THEN
         NFR = N0
         I1 = 1
      ELSE
         NFR = NFT+1
         I1 = -1
      ENDIF
*
      DO 4 I = 1,NFT
         IK = NFR+I1*I
         IT = ITFREE(IK)
         K = NBT1+IK
*
         DO 3 J = 1,3
            JDIAG(K,J) = JDIAG(IT,J)
            ARR(K,J) = ARR(IT,J)
    3    CONTINUE
*
    4 CONTINUE
*
      GOTO (22,23,24),ISW
*
      END
