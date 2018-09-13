************************************************************************
*                                                                      *
      SUBROUTINE SPRATE (M)
*                                                                      *
*   This  subroutine  prepares  the  information to be transfered to   *
*   GENSUM for numerical evaluation.                                   *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*6 NAME
      CHARACTER*5 NME
      LOGICAL SUM6J,T6J,JT,JS,SUMVAR,CUT
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION JTEM4(MTRIAD,M6J),JTEM5(MTRIAD,M6J),JTEM6(MTRIAD),
     :   NSUM6J(M6J),J6SUM(M6J)
      DIMENSION SUM6J(M6J),T6J(M6J),JT(MTRIAD),JS(MTRIAD),
     :   INVER(MANGM),JNSUM(MTRIAD),JINV(MTRIAD),N6JN(M6J),IN6J(M6J),
     :   JSUMT(M6J,6)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),JW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CUTDIG/CUT
     :      /DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
     :      /SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :              JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :              K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :              JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
*
*   Test that array dimensions have not been exceeded.
*
      IF (MP .GT. MANGM) THEN
         NMX = MANGM
         NPX = MP
         NAME = 'MANGM '
         NME  = 'MP   '
      ELSEIF (JWC .GT. M6J) THEN
         NMX = M6J
         NPX = JWC
         NAME = 'M6J   '
         NME  = 'JWC  '
      ELSEIF (J6C .GT. M3MNGM) THEN
         NMX = M3MNGM
         NPX = J6C
         NAME = 'M3MNGM'
         NME  = 'J6C  '
      ELSEIF (J7C .GT. M3MNGM) THEN
         NMX = M3MNGM
         NPX = J7C
         NAME = 'M3MNGM'
         NME  = 'J7C  '
      ELSEIF (J8C .GT. M3MNGM) THEN
         NMX = M3MNGM
         NPX = J8C
         NAME = 'M3MNGM'
         NME  = 'J8C  '
      ELSE
         IF (J9C .LE. MANGMP) GOTO 54
         NMX = MANGMP
         NPX = J9C
         NAME = 'MANGMP'
         NME  = 'J9C  '
      ENDIF
*
   60 WRITE (*,300) NAME,NME,NPX,NMX
      STOP
*
*   Determination of effective summation variables and their
*   relationships with 6j coefficients.
*
   54 DO 2 I = 1,JWC
         INV6J(I) = 0
         SUM6J(I) = .FALSE.
    2 CONTINUE
*
      NSUM = 0
      NLSUM = 0
      IF (MP .EQ. M) RETURN
      M1 = M+1
*
      DO 1 I = M1,MP
         IF (SUMVAR(I)) THEN
            NSUM = NSUM+1
            JSUM6(NSUM) = 0
            INVER(I) = NSUM
         ENDIF
    1 CONTINUE
*
      IF (NSUM .EQ. 0) RETURN
*
      IF (NSUM .GT. MTRIAD) THEN
         NMX = MTRIAD
         NPX = NSUM
         NAME = 'MTRIAD'
         NME  = 'NSUM '
         GOTO 60
      ENDIF
*
      KT = 0
*
      DO 4 I = 1,JWC
         DO 5 J = 1,6
            IK = JW(J,I)
            IF (.NOT. SUMVAR(IK)) GOTO 5
*
            IF (.NOT. SUM6J(I)) THEN
               SUM6J(I) = .TRUE.
               KT = KT+1
               J6SUM(KT) = 0
               NSUM6J(KT) = I
               INV6J(I) = KT
            ENDIF
*
            ISK = INVER(IK)
            I2 = JSUM6(ISK)+1
            JSUM6(ISK) = I2
            JSUM4(ISK,I2) = J
            JSUM5(ISK,I2) = KT
            I3 = J6SUM(KT)+1
            J6SUM(KT) = I3
            JSUMT(KT,I3) = ISK
    5    CONTINUE
    4 CONTINUE
*
      CALL VAR (J6,J6P,J6C,J6CP,J6CC,SUMVAR,MP,M,INVER)
      CALL VAR (J7,J7P,J7C,J7CP,J7CC,SUMVAR,MP,M,INVER)
      CALL VAR (J8,J8P,J8C,J8CP,J8CC,SUMVAR,MP,M,INVER)
      CALL VAR (J9,J9P,J9C,J9CP,J9CC,SUMVAR,MP,M,INVER)
*
      IF (.NOT. CUT) THEN
         NLSUM = 1
         NBJ(1) = NSUM
         NB6J(1) = KT
         K6CP(1) = J6CP
         K7CP(1) = J7CP
         K8CP(1) = J8CP
         K9CP(1) = J9CP
*
         DO 21 I = 1,KT
            I1 = NSUM6J(I)
            DO 22 J = 1,6
               JWORD(J,I) = JW(J,I1)
   22       CONTINUE
   21    CONTINUE
*
         DO 80 I = 1,NSUM
            ISU = JSUM6(I)
            DO 81 J = 1,ISU
               I1 = JSUM5(I,J)
               J1 = JSUM4(I,J)
               JWORD(J1,I1) = MP+I
   81       CONTINUE
   80    CONTINUE
*
         RETURN
      ENDIF
*
*   Separation of variables and sums in case a cut was detected.
*
      K6C = 0
      K7C = 0
      K8C = 0
      K9C = 0
      NJ = 0
      N6J = 0
*
      DO 9 I = 1,KT
         T6J(I) = .FALSE.
    9 CONTINUE
*
      DO 7 I = 1,NSUM
         JT(I) = .FALSE.
         JS(I) = .FALSE.
    7 CONTINUE
*
      J = 1
*
   10 NJ = NJ+1
      JNSUM(NJ) = J
      JINV(J) = NJ
      JT(J) = .TRUE.
   18 JS(J) = .TRUE.
      JS6 = JSUM6(J)
*
      DO 11 I = 1,JS6
         I6J = JSUM5(J,I)
*
         IF (.NOT. T6J(I6J)) THEN
            T6J(I6J) = .TRUE.
            N6J = N6J+1
            N6JN(N6J) = NSUM6J(I6J)
            IN6J(I6J) = N6J
         ENDIF
*
         J6J = J6SUM(I6J)
*
         DO 12 K = 1,J6J
            JK = JSUMT(I6J,K)
            IF (.NOT. JT(JK)) THEN
               NJ = NJ+1
               JNSUM(NJ) = JK
               JINV(JK) = NJ
               JT(JK) = .TRUE.
            ENDIF
   12    CONTINUE
*
   11 CONTINUE
*
      DO 13 JJ = 1,NSUM
         J = JJ
         IF ((.NOT. JS(JJ)) .AND. JT(JJ)) GOTO 18
   13 CONTINUE
*
      NLSUM = NLSUM+1
*
      IF (NLSUM .GT. MSUM) THEN
         NMX = MSUM
         NPX = NLSUM
         NAME = 'MSUM  '
         NME = 'NLSUM'
         GOTO 60
      ENDIF
      NBJ(NLSUM) = NJ
      NB6J(NLSUM) = N6J
*
      IF (J6CP .NE. 0) CALL CHVAR (J6P,J6CP,K6C,JT,JINV,NSUM)
      K6CP(NLSUM) = K6C
      IF (J7CP .NE. 0) CALL CHVAR (J7P,J7CP,K7C,JT,JINV,NSUM)
      K7CP(NLSUM) = K7C
      IF (J8CP .NE. 0) CALL CHVAR (J8P,J8CP,K8C,JT,JINV,NSUM)
      K8CP(NLSUM) = K8C
      IF (J9CP .NE. 0) CALL CHVAR (J9P,J9CP,K9C,JT,JINV,NSUM)
      K9CP(NLSUM) = K9C
*
      IF (NJ .NE. NSUM) THEN
         DO 16 JJ = 1,NSUM
            J = JJ
            IF (.NOT. JT(JJ)) GOTO 10
   16    CONTINUE
      ENDIF
*
      DO 26 I = 1,KT
         I1 = N6JN(I)
         DO 27 J = 1,6
            JWORD(J,I) = JW(J,I1)
   27    CONTINUE
   26 CONTINUE
*
      DO 28 I = 1,NSUM
         IK = JNSUM(I)
         I2 = JSUM6(IK)
         JTEM6(I) = I2
         DO 29 J = 1,I2
            JTEM4(I,J) = JSUM4(IK,J)
            K = JSUM5(IK,J)
            JTEM5(I,J) = IN6J(K)
   29    CONTINUE
   28 CONTINUE
*
      DO 40 I = 1,NSUM
         I2 = JTEM6(I)
         JSUM6(I) = I2
         DO 41 J = 1,I2
            I1 = JTEM5(I,J)
            J1 = JTEM4(I,J)
            JSUM4(I,J) = J1
            JSUM5(I,J) = I1
            JWORD(J1,I1) = I+MP
   41    CONTINUE
   40 CONTINUE
*
      RETURN
*
  300 FORMAT (' Dimension error for ',A6
     :       /2X,A5,' = ',I5,' is out of allowed range',I4)
*
      END
