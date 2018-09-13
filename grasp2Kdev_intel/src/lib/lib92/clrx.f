************************************************************************
*                                                                      *
      FUNCTION CLRX (KAPPAA,K,KAPPAB)
*                                                                      *
*   The value of CLRX is the 3-j symbol:                               *
*                                                                      *
*                    ( JA        K        JB  )                        *
*                    ( 1/2       0       -1/2 )                        *
*                                                                      *
*   The  K'S are kappa angular quantum numbers. The formula is taken   *
*   from D M Brink and G R Satchler, <Angular Momentum>, second edi-   *
*   tion (Oxford: Clarendon press, 1968), p 138.   The logarithms of   *
*   the first  MFACT  factorials must be available in  COMMON/FACTS/   *
*   for this program to function correctly. Note that  N!  is stored   *
*   in FACT(N+1)                                                       *
*                                                                      *
*   No subroutines called.                                             *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
*
      PARAMETER (MFACT = 500)
*
      COMMON/FACTS/GAM(MFACT)
*
*   Determine the absolute values of the kappas
*
      KA = ABS (KAPPAA)
      KB = ABS (KAPPAB)
*
*   Perform the triangularity check
*
      IF ((ABS (KA-KB) .LE. K) .AND. (KA+KB-1 .GE. K)) THEN
*
*   Triangularity satisfied; compute the 3j coefficient
*
*   Begin with the logarithm of the square of the leading term
*
         EXPTRM = -LOG (DBLE (KA*KB))
*
*   Compute the logarithm of the square root of the leading term
*   and the factorial part that doesn't depend on the parity of
*   KA+KB+K (the delta factor)
*
         KAPKB = KA+KB
         KABKP = KAPKB+K
         KAMKB = KA-KB
         KBMKA = KB-KA
         EXPTRM = 0.5D 00
     :           *(EXPTRM+GAM(KAPKB-K  )+GAM(KAMKB+K+1)
     :                   +GAM(KBMKA+K+1)-GAM(KABKP  +1) )
*
*   The remainder depends on the parity of KA+KB+K
*
         IF (MOD (KABKP,2) .EQ. 0) THEN
*
*   Computation for even parity case
*
*   Include the phase factor: a minus sign if necessary
*
            IF (MOD (3*KABKP/2,2) .EQ. 0) THEN
               CLRX =  1.0D 00
            ELSE
               CLRX = -1.0D 00
            ENDIF
*
*   Include the contribution from the factorials
*
            EXPTRM = EXPTRM+GAM((KABKP  +2)/2)-GAM((KAPKB-K  )/2)
     :                     -GAM((KAMKB+K+2)/2)-GAM((KBMKA+K+2)/2)
*
         ELSE
*
*   Computation for odd parity case
*
*   Include the phase factor: a minus sign if necessary
*
            IF (MOD ((3*KABKP-1)/2,2) .EQ. 0) THEN
               CLRX =  1.0D 00
            ELSE
               CLRX = -1.0D 00
            ENDIF
*
*   Include the contribution from the factorials
*
            EXPTRM = EXPTRM+GAM((KABKP  +1)/2)-GAM((KAPKB-K+1)/2)
     :                     -GAM((KAMKB+K+1)/2)-GAM((KBMKA+K+1)/2)
*
         ENDIF
*
*   Final assembly
*
         CLRX = CLRX*EXP (EXPTRM)
*
      ELSE
*
*   Triangularity violated; set the coefficient to zero
*
         CLRX = 0.0D 00
*
      ENDIF
*
      RETURN
*
      END
