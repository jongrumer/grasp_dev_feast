************************************************************************
*                                                                      *
      SUBROUTINE SETTAB (FAIL)
*                                                                      *
*   Builds up the unstructured graph. Sets the array J23, containing   *
*   the  two lists of original triads J2 and J3, and the correspond-   *
*   ing arrows  on the  angular  momenta lines. Also establishes the   *
*   numerical and phase factors  connecting  recoupling  coefficient   *
*   and graphs, according to Yutsis, Levinson, and Vanagas. For this   *
*   purpose determines the total J.                                    *
*                                                                      *
*                                           Last update: 16 Ocy 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      INTEGER ARROW
      LOGICAL FAIL,TABS,FREE,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
     :      /BUILD/IAL(M4TRD),IF1,IF2,NODE
     :      /ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
*
      DATA NAME/'SETTAB'/
*
      DO 502 I = 1,M2TRD
         DO 501 J = 1,3
            J23(I,J) = 0
  501    CONTINUE
  502 CONTINUE
*
      IPR = N-1
      NBTR = IPR+IPR
*
      DO 4 I = 1,IPR
         DO 5 J = 1,2
            J23(I,J) = J2(I,J)
            ARROW(I,J) = 1
    5    CONTINUE
         TABS(I) = .FALSE.
         J23(I,3) = J2(I,3)
         ARROW(I,3) = -1
    4 CONTINUE
*
      IPR1 = IPR+1
*
      DO 7 I = IPR1,NBTR
         II = I-IPR
         DO 6 J = 1,2
            J23(I,J) = J3(II,J)
            ARROW(I,J) = -1
    6    CONTINUE
         TABS(I) = .FALSE.
         J23(I,3) = J3(II,3)
         ARROW(I,3) = 1
    7 CONTINUE
*
      DO 11 J = 1,NBTR
         J8(J) = J23(J,1)
   11 CONTINUE
*
      J8C = NBTR+IPR
      NB1 = NBTR+1
*
      DO 12 J = NB1,J8C
         I = J-IPR
         J8(J) = J23(I,3)
   12 CONTINUE
*
      J6C = NBTR
*
      DO 13 J = 1,J6C
         J6(J) = J23(J,3)
   13 CONTINUE
*
      DO 10 I = 1,M
         SUMVAR(I) = .FALSE.
         IAL(I) = 1
   10 CONTINUE
*
      DO 9 I = 1,NBTR
         DO 8 J = 1,3
            JI = J23(I,J)
            K = IAL(JI)
            LINE(JI,K) = I
            LCOL(JI,K) = J
            IAL(JI) = K+1
    8    CONTINUE
    9 CONTINUE
*
      IT = 0
*
      DO 18 I = 1,NBTR
*
         JT = J23(I,3)
*
         IF (IAL(JT) .EQ. 3) THEN
*
            CALL OTHERJ (I,JT,L,LC,K)
            IF (LC .EQ. 3) GOTO 19
*
         ELSE
*
            IF (IT .EQ. 1) THEN
               CALL DELTA (JT1,JT,FAIL)
               IF (FAIL) GOTO 20
               K = LINE(JT,1)
               KC = LCOL(JT,1)
               LINE(JT1,2) = K
               LCOL(JT1,2) = KC
               LINE(JT,2) = LINE(JT1,1)
               LCOL(JT,2) = LCOL(JT1,1)
               J23(K,KC) = JT1
               IAL(JT) = 1
               GOTO 19
            ENDIF
*
            JT1 = JT
            IT = 1
*
         ENDIF
*
   18 CONTINUE
*
   19 J9(J9C+1) = JT
      J9C = J9C+2
      J9(J9C) = JT
*
   20 CALL PRINTJ (NAME,4)
*
      RETURN
      END
