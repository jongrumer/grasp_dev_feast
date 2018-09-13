************************************************************************
*                                                                      *
      SUBROUTINE TRDEL (JJ1,JJ2,JJ3,NBN,FAIL)
*                                                                      *
*   Test for triangular delta. If not satisfied FAIL = .TRUE. .        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL FAIL,SUMVAR,CUT,FREE
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :            J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),
     :            MP
     :      /CUTDIG/CUT
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      IF (SUMVAR(JJ1) .OR. SUMVAR(JJ2) .OR. SUMVAR(JJ3)) RETURN
      IF (NBN .GT. 4) CUT = .TRUE.
      IF ((.NOT. FREE(JJ1)) .AND. (.NOT. FREE(JJ2)) .AND.
     :    (.NOT.FREE(JJ3))) THEN
         I1 = J1(JJ1)
         I2 = J1(JJ2)
         I3 = J1(JJ3)
         IF ((I1 .LT. (ABS (I2-I3)+1)) .OR. (I1 .GT. (I2+I3-1)))
     :      FAIL = .TRUE.
      ENDIF
*
      RETURN
      END
