************************************************************************
*                                                                      *
      SUBROUTINE KNJ (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,KDW,
     :            JDDEL,LDDEL,DSUMVR,MDP,
     :            JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     :            KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,INVD6J)
*                                                                      *
*   This routine stores data for future calls to GENSUM.               *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      PARAMETER (MANGM = 60, MTRIAD = 12)
      PARAMETER (M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3))
      PARAMETER (M6J = 20, MSUM = 10)
*
      LOGICAL SUMVAR,DSUMVR
*
      DIMENSION JD6(M3MNGM),JD7(M3MNGM),JD8(M3MNGM),
     :          JD9(MANGMP),KDW(6,M6J),LDDEL(M6J,2),DSUMVR(MANGM)
      DIMENSION JD6P(MANGMP),JD7P(MANGMP),JD8P(MANGMP),JD9P(MANGMP),
     :          JDWORD(6,M6J),NDBJ(MSUM),
     :          NDB6J(MSUM),KD6CP(MSUM),KD7CP(MSUM),KD8CP(MSUM),
     :          KD9CP(MSUM),JDSUM6(MTRIAD),JDSUM4(MTRIAD,M6J),
     :          JDSUM5(MTRIAD,M6J),INVD6J(M6J)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :       J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),MP
     :      /SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :       JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :       K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :       JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
*
      JD6C = J6C
      JD7C = J7C
      JD8C = J8C
      JD9C = J9C
      JDWC = JWC
      JDDEL = JDEL
      MDP = MP
      NDLSUM = NLSUM
      IF (J6C .NE. 0) THEN
         DO 1 I = 1,J6C
            JD6(I) = J6(I)
  1      CONTINUE
      ENDIF
      IF (J7C .NE. 0) THEN
         DO 2 I = 1,J7C
            JD7(I) = J7(I)
  2      CONTINUE
      ENDIF
      IF (J8C .NE. 0) THEN
         DO 3 I = 1,J8C
            JD8(I) = J8(I)
  3      CONTINUE
      ENDIF
      IF (J9C .NE. 0) THEN
         DO 4 I = 1,J9C
            JD9(I) = J9(I)
  4      CONTINUE
      ENDIF
      IF (JWC .NE. 0) THEN
         DO 5 I = 1,6
            DO 5 J = 1,JWC
               KDW(I,J) = KW(I,J)
  5      CONTINUE
         DO 6 I = 1,JWC
            INVD6J(I) = INV6J(I)
  6      CONTINUE
      ENDIF
      IF (JDEL .NE. 0) THEN
         DO 7 I = 1,2
            DO 7 J = 1,JDEL
               LDDEL(J,I) = LDEL(J,I)
  7      CONTINUE
      ENDIF
      IF (MP .NE. 0) THEN
         DO 8 I = 1,MP
            DSUMVR(I) = SUMVAR(I)
  8      CONTINUE
      ENDIF
      IF (NLSUM .NE. 0) THEN
         DO 9 I = 1,NLSUM
            NDBJ(I) = NBJ(I)
            NDB6J(I) = NB6J(I)
            KD6CP(I) = K6CP(I)
            KD7CP(I) = K7CP(I)
            KD8CP(I) = K8CP(I)
            KD9CP(I) = K9CP(I)
  9      CONTINUE
      ENDIF
      DO 10 I = 1,MANGMP
         JD6P(I) = J6P(I)
         JD7P(I) = J7P(I)
         JD8P(I) = J8P(I)
         JD9P(I) = J9P(I)
  10  CONTINUE
      DO 11 I = 1,MTRIAD
         JDSUM6(I) = JSUM6(I)
         DO 11 J = 1,M6J
            JDSUM4(I,J) = JSUM4(I,J)
            JDSUM5(I,J) = JSUM5(I,J)
  11  CONTINUE
      DO 12 I = 1,6
         DO 12 J = 1,M6J
            JDWORD(I,J) = JWORD(I,J)
  12  CONTINUE
*
      RETURN
      END
