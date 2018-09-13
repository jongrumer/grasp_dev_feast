************************************************************************
*                                                                      *
      SUBROUTINE SIGMA_1 (IPAR,I1,I2,APART)
*                                                                      *
*   This  routine computes                                             *
*   [-kappa_a || sigma^(1) ||kappa_b]  if IPAR = 1                     *
*   <-kappa_a || sigma^(1) ||kappa_b>  if IPAR = 2                     *
*                                                                      *
*                                                                      *
*   Written by G. Gaigalas and E. Gaidamauskas                         *
*                                         Last revision: 16 Doc 2009   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
*
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
*
*   Set KAP1 and KAP2
*
      KAP1 = NAK(I1)
      KAP2 = NAK(I2)
*
*   Determine the l quantum numbers
*
      IF (KAP1 .GT. 0) THEN
         L1 =  KAP1
      ELSE
         L1 = -KAP1-1
      ENDIF
*
      IF (KAP2 .GT. 0) THEN
         L2 =  KAP2
      ELSE
         L2 = -KAP2-1
      ENDIF
*
*   Determine the j quantum numbers and l_1 tilda
*
      J1 = IABS(NAK(I1))*2-1
      J2 = IABS(NAK(I2))*2-1
      L1_TILDA = J1 - L1
      IF (L1_TILDA .EQ. L2) THEN
*
*   Determine the Racah W coefficients.
*
         CALL DRACAH (1,J1,1,J2,2*L1_TILDA,2,RAC)
         IF (MOD(1+J1+1+J2,4).EQ.0) RAC=-RAC
         APART =  RAC*DSQRT(6.0D 00*(J2+1)) 
      ELSE
         APART = 0.0D 00
      END IF
*
*   Determine the phase factori
*
      IF(MOD(2*L1_TILDA+1+J1+2,4).EQ.0) APART=-APART
      IF(IPAR .EQ. 2) THEN
         APART = DSQRT(DBLE(J1+1))*APART
      END IF
      RETURN
      END
