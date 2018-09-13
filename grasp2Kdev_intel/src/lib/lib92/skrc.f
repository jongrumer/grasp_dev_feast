************************************************************************
*                                                                      *
      SUBROUTINE SKRC (IS,KAPS,KS,KD1,KD2,KE1,KE2)
*                                                                      *
*   Determines the range of the tensor rank k for Coulomb integral.    *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      DIMENSION IS(4),KAPS(4),KS(4)
*
      KD2 = 0
      KE2 = 0
*
*   Direct terms --- KD1 = minimum k , KD2 = number of terms
*
      ISD1 = 1
      IF (KAPS(1)*KAPS(3) .LT. 0) ISD1 = -1
      ISD2 = 1
      IF (KAPS(2)*KAPS(4) .LT. 0) ISD2 = -1
      KD1A = ABS (KS(1)-KS(3))
      IF (ISD1 .LT. 0) KD1A = KD1A+2
      KD1B = ABS (KS(2)-KS(4))
      IF (ISD2 .LT. 0) KD1B = KD1B+2
      IF (MOD ((KD1A-KD1B)/2,2) .NE. 0) GOTO 1
      KD2A = KS(1)+KS(3)-2
      IF (ISD1 .GT. 0) KD2A = KD2A-2
      KD2B = KS(2)+KS(4)-2
      IF (ISD2 .GT. 0) KD2B = KD2B-2
      KD1 = MAX (KD1A,KD1B)/2
      KD2 = MIN (KD2A,KD2B)/2
      KD2 = (KD2-KD1)/2+1
*
*   Exchange terms --- KE1 = minimum k , KE2 = number of terms
*
    1 CONTINUE
      IF ((IS(1) .EQ. IS(2)) .OR. (IS(3) .EQ. IS(4))) RETURN
      ISE1 = 1
      IF (KAPS(1)*KAPS(4) .LT. 0) ISE1 = -1
      ISE2 = 1
      IF (KAPS(2)*KAPS(3) .LT. 0) ISE2 = -1
      KE1A = ABS (KS(1)-KS(4))
      IF (ISE1 .LT. 0) KE1A = KE1A+2
      KE1B = ABS (KS(2)-KS(3))
      IF (ISE2 .LT. 0) KE1B = KE1B+2
      IF (MOD ((KE1A-KE1B)/2,2) .NE. 0) RETURN
      KE2A = KS(1)+KS(4)-2
      IF (ISE1 .GT. 0) KE2A = KE2A-2
      KE2B = KS(2)+KS(3)-2
      IF (ISE2 .GT. 0) KE2B = KE2B-2
      KE1 = MAX (KE1A,KE1B)/2
      KE2 = MIN (KE2A,KE2B)/2
      KE2 = (KE2-KE1)/2+1
*
      RETURN
      END
