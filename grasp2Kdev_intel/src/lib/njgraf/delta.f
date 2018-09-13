************************************************************************
*                                                                      *
      SUBROUTINE DELTA (JA,JB,FAIL)
*                                                                      *
*   Test for delta(JA,JB). If they are summation variables, the sec-   *
*   ond  is  changed  into  the first everywhere. if they are fixed,   *
*   their value is checked, and fail put to .TRUE. if they differ.     *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL FAIL,CUT,SUMVAR,FREE
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CUTDIG/CUT
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DIM/J6CC,J7CC,J8CC,J9CC,JWCC,JDELC
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      IF (IBUG3 .EQ. 1) WRITE (99,300) JA,SUMVAR(JA),JB,SUMVAR(JB)
      IF (SUMVAR(JA) .AND. SUMVAR(JB)) GOTO 2
      IF (FREE(JA) .OR. FREE(JB)) THEN
         JDEL = JDEL+1
         LDEL(JDEL,1) = JA
         LDEL(JDEL,2) = JB
         SUMVAR(JA) = .FALSE.
         SUMVAR(JB) = .FALSE.
         RETURN
      ENDIF
*
      IF (J1(JA) .NE. J1(JB)) FAIL = .TRUE.
      CUT = .TRUE.
      RETURN
*
    2 IF (J6C .NE. J6CC) THEN
         J61 = J6CC+1
*
         DO 3 I = J61,J6C
            IF (J6(I) .EQ. JB) J6(I) = JA
    3    CONTINUE
*
      ENDIF
*
      IF (J7C .NE. J7CC) THEN
         J71 = J7CC+1
*
         DO 5 I = J71,J7C
            IF (J7(I) .EQ. JB) J7(I) = JA
    5    CONTINUE
      ENDIF
*
      IF (J8C .NE. J8CC) THEN
         J81 = J8CC+1
*
         DO 7 I = J81,J8C
            IF (J8(I) .EQ. JB) J8(I) = JA
    7    CONTINUE
      ENDIF
*
      IF (J9C .NE. J9CC) THEN
         J91 = J9CC+1
*
         DO 9 I = J91,J9C
            IF (J9(I) .EQ. JB) J9(I) = JA
    9    CONTINUE
      ENDIF
*
      IF (JWC .NE. JWCC) THEN
         JW1 = JWCC+1
*
         DO 14 I = JW1,JWC
            DO 13 J = 1,6
               IF (KW(J,I) .EQ. JB) KW(J,I) = JA
   13       CONTINUE
   14    CONTINUE
      ENDIF
*
      IF (JDEL .NE. JDELC) THEN
         JDEL1 = JDELC+1
*
         DO 17 I = JDEL1,JDEL
            DO 16 J = 1,2
               IF (LDEL(I,J) .EQ. JB) LDEL(I,J) = JA
   16       CONTINUE
   17    CONTINUE
*
         SUMVAR(JB) = .FALSE.
      ENDIF
*
      RETURN
*
  300 FORMAT (/'From DELTA: JA = ',I2,L2,5X,'JB = ',I2,L2)
*
      END
