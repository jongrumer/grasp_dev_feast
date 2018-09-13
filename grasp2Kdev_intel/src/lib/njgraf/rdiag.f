************************************************************************
*                                                                      *
      SUBROUTINE RDIAG (I,J,IK1,IK2,ICHAN,MAT,JSUM,J12)
*                                                                      *
*   Called by  GENSUM to establish the range of values of the summa-   *
*   tion variables.  This routine replaces an extended range do loop   *
*   in GENSUM, to conform with the FORTRAN 77 standard.                *
*                                                                      *
*                                           Last update: 02 Sep 1987   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      PARAMETER (MTRIAD = 12,M6J = 20)
*
      DIMENSION MAT(MTRIAD,MTRIAD),JMNP(5),JMXP(5)
      DIMENSION J12(4,MTRIAD,MTRIAD),JSUM(2,M6J)
*
      JMIN1 = 0
      JMAX1 = 1000
      K = MAT(IK1,IK2)
*
      DO 1 L1 = 1,K
*
         L3 = MAT(J,J)
         JJ1 = JSUM(1,J)
         JND = J12(L1,IK2,IK1)
         JMIN = 1000
         JMAX = 0
         JMNP(L1) = 0
         JMXP(L1) = 1000
*
         DO 2 L2 = 1,L3
*
            JMN = IABS(JND-JJ1)
            JMX = JND+JJ1
            JMIN = MIN (JMN,JMIN)
            JMAX = MAX (JMX,JMAX)
            JMNP(L1) = MAX (JMN,JMNP(L1))
            JMXP(L1) = MIN (JMX,JMXP(L1))
            JJ1 = JJ1+2
*
    2    CONTINUE
*
         JMIN1 = MAX (JMIN1,JMIN)
         JMAX1 = MIN (JMAX1,JMAX)
*
    1 CONTINUE
*
      IF (MAT(I,I) .EQ. 0) THEN
         JSUM(1,I) = JMIN1
         JSUM(2,I) = JMAX1
         MAT(I,I) = (JMAX1-JMIN1)/2+1
         ICHAN = ICHAN+1
         GOTO 3
      ENDIF
*
      IF (JSUM(1,I) .LT. JMIN1) THEN
         JSUM(1,I) = JMIN1
         ICHAN = ICHAN+1
      ENDIF
*
      IF (JSUM(2,I) .GT. JMAX1) THEN
         JSUM(2,I) = JMAX1
         ICHAN = ICHAN+1
      ENDIF
*
    3 K1 = 0
*
      DO 4 L1 = 1,K
         IF ((JMNP(L1) .LE. JSUM(1,I)) .AND.
     :       (JMXP(L1) .GE. JSUM(2,I))) GOTO 4
         K1 = K1+1
         J12(K1,IK2,IK1) = J12(L1,IK2,IK1)
    4 CONTINUE
*
      IF (K1 .NE. K) THEN
         MAT(IK1,IK2) = K1
         ICHAN = ICHAN+1
      ENDIF
*
      MAT(IK2,IK1) = J12(1,IK2,IK1)
*
      RETURN
      END
