************************************************************************
*                                                                      *
      SUBROUTINE VAR (JN,JNS,JNC,JNSC,JBC,SUMVAR,MP,M,INVER)
*                                                                      *
*   Test  for  variable  character and put in JNS if yes, and JN now   *
*   contains 0.                                                        *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL SUMVAR(MP)
*
      PARAMETER (
     :   MANGM = 60,M3MNGM = 3*MANGM,MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12,M2TRD = 2*MTRIAD,M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10)
*
      DIMENSION JN(JNC),JNS(MANGMP),INVER(MP)
*
      JNSC = 0
      IF (JBC .NE. JNC) THEN
         JBBC = JBC+1
*
         DO 1 I = JBBC,JNC
            I1 = JN(I)
            IF (SUMVAR(I1)) THEN
               JNSC = JNSC+1
               IF (JNSC .GT. MANGMP) THEN
                  WRITE (*,300) JNSC,MANGMP
                  STOP
               ENDIF
               J = INVER(I1)
               JNS(JNSC) = J
               JN(I) = M
            ENDIF
    1    CONTINUE
      ENDIF
*
      RETURN
*
  300 FORMAT (' Dimension error in VAR. JNSC = ',I5,' MANGMP = ',I5)
*
      END
