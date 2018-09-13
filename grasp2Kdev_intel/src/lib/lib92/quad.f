************************************************************************
*                                                                      *
      SUBROUTINE QUAD (RESULT)
*                                                                      *
*   The argument result is an approximation  to the integral of F(R)   *
*   from  zero  to  infinity,  where the values of RP(I)*F(R(I)) are   *
*   tabulated in the array  TA(I). The integral in the interval zero   *
*   to R(J) is computed by use of an analytical fit                    *
*                                                                      *
*                                SIGMA                                 *
*                      F(R) = A R                                      *
*                                                                      *
*   A five-point  Closed  Newton-Cotes  formula (cf. F B Hildebrand,   *
*   Introduction to Numerical Analysis, second edition, McGraw-Hill,   *
*   New York, 1974, p 93)  is  used  to  compute the integral in the   *
*   interval  R(J:MTP).  The  contribution  from  the  tail  of  the   *
*   function beyond the last  tabular  point  (MTP) is assumed to be   *
*   negligible. The method uses  MTP+3  tabulation points. Array  TA   *
*   should therefore be dimensioned to at least  N+4 .                 *
*                                                                      *
*   No subroutines called.                                             *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      COMMON/DEF0/TENMAX,EXPMAX,EXPMIN,PRECIS
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NCC/C1,C2,C3,C4
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Find first values that will permit computation of exponent
*
      MTPM1 = MTP-1
      DO 3 I = 2,MTPM1
*
         TAI = TA(I)
         IF (ABS (TAI) .GT. 0.D0) THEN
*
            IP1   = I+1
            TAIP1 = TA(IP1)
            QUOTT = TAIP1/TAI
*
            IF (QUOTT .GT. 0.D0) THEN
*
*   Exponent from fit
*
               FRIP1 = TAIP1/RP(IP1)
               FRI   = TAI  /RP(I  )
               RATIO = FRIP1/FRI
               RIP1  = R (IP1)
               RI    = R (I  )
               SIGMA = LOG (RATIO)/LOG (RIP1/RI)
*
*   Analytical integration and error estimate for interval r(1:i)
*
               FRI    = RI*FRI
               RESULT = FRI/(SIGMA+1.D0)
*
*   Set the tail to zero
*
               DO 1 LOC = 1,3
                  TA(MTP+LOC) = 0.D0
    1          CONTINUE
*
*   Newton-Cotes quadature for the remainder
*
               RESULT = RESULT+C1*TAI
               DO 2 LOC = IP1,MTP,4
                  RESULT = RESULT+C2*(TA(LOC  )+TA(LOC+2))
     :                           +C3* TA(LOC+1)
     :                           +C4* TA(LOC+3)
    2          CONTINUE
               IF (MOD (MTP-I,4) .EQ. 0) RESULT = RESULT-C1*TA(MTP)
*
*   Test of result's accuracy; `decomment' to activate
*
*              ESTDER = 10.0D 00*RI*FRI
*              RATIO = ABS (ESTDER/RESULT)
*              IF (RATIO .GT. ACCY) PRINT (*,300) RATIO
*
               GOTO 4
*
            ENDIF
*
         ENDIF
*
    3 CONTINUE
*
*   No value which will permit computation of exponent
*
      RESULT = 0.D0
*
    4 RETURN
*
  300 FORMAT (/'QUAD: Estimated accuracy is ',1PD10.3,
     :        /' Decrease RNT or improve input data conditioning to'
     :        ,' ameliorate.'/)
*
      END
