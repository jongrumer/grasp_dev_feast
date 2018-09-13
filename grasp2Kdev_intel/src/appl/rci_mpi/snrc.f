************************************************************************
*                                                                      *
      SUBROUTINE SNRC (IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
*                                                                      *
*   Determines the range of tensor rank NU for direct/exchange terms,  *
*   and classifies the types of radial integral.                       *
*                                                                      *
*   Input variables:                                                   *
*                                                                      *
*      IS      : Orbital labels                                        *
*      KAPS    : Values of 2*kappa                                     *
*      KS      : Values of 2*J+1                                       *
*                                                                      *
*   Outputs:                                                           *
*                                                                      *
*      ND1/NE1 : Lowest NU value for direct/exchange types             *
*      ND2/NE2 : Corresponding number of contributing NU values: NU    *
*                = ND1, ND1+2 ,..., ND1+2*(ND2-1) etc                  *
*      IBRD    : Classify types of  radial  integrals  contributing;   *
*      IBRE      negative value implies null contribution              *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      DIMENSION IS(4),KAPS(4),KS(4)
*
      ND2 = 0
      NE2 = 0
*
*   2.0  Form limits for direct terms
*
      IAC = 1
      IF ((KAPS(1)*KAPS(3)) .LT. 0) IAC = -1
      IAD = 1
      IF ((KAPS(2)*KAPS(4)) .LT. 0) IAD = -1
      ND1 = ABS (KS(1)-KS(3))/2-1
      IF (IAC .EQ. -1) ND1 = ND1+1
      IF (ND1 .EQ. -1) ND1 = 1
      ND1A = ABS (KS(2)-KS(4))/2-1
      IF (IAD .EQ. -1) ND1A = ND1A+1
      IF (ND1A .EQ. -1) ND1A = 1
      IF (MOD (ND1-ND1A,2) .EQ. 0) GOTO 1
      IBRD = -1
      GOTO 2
    1 ND2 = ABS (KS(1)+KS(3))/2
      IF (IAC .EQ. -1) ND2 = ND2+1
      ND2A = ABS (KS(2)+KS(4))/2
      IF (IAD .EQ. -1) ND2A = ND2A+1
      ND1 = MAX (ND1,ND1A)
      ND2 = MIN (ND2,ND2A)
      ND2 = (ND2-ND1)/2+1
*
*   2.1  Identify type of radial integrals
*
      IBRD = 1
      IF (IS(1) .EQ. IS(3) .AND. IS(2) .NE. IS(4)) IBRD = 2
      IF (IS(1) .NE. IS(3) .AND. IS(2) .EQ. IS(4)) IBRD = 2
      IF (IS(1) .EQ. IS(3) .AND. IS(2) .EQ. IS(4)) IBRD = 3
*
*   3.0  Form limits for exchange terms
*
    2 IF ((IS(1) .NE. IS(2)) .AND. (IS(3) .NE. IS(4))) GOTO 3
      IBRE = -1
      RETURN
    3 CONTINUE
      IAC = 1
      IF ((KAPS(1)*KAPS(4)) .LT. 0) IAC = -1
      IAD = 1
      IF ((KAPS(2)*KAPS(3)) .LT. 0) IAD = -1
      NE1 = IABS (KS(1)-KS(4))/2-1
      IF (IAC .EQ. -1) NE1 = NE1+1
      IF (NE1 .EQ. -1) NE1 = 1
      NE1A = ABS (KS(2)-KS(3))/2-1
      IF (IAD .EQ. -1) NE1A = NE1A+1
      IF (NE1A .EQ. -1) NE1A = 1
      IF (MOD (NE1-NE1A,2) .EQ. 0) GOTO 4
      IBRE = -1
      RETURN
*
    4 NE2 = ABS (KS(1)+KS(4))/2
      IF (IAC .EQ. -1) NE2 = NE2+1
      NE2A = ABS (KS(2)+KS(3))/2
      IF (IAD .EQ. -1) NE2A = NE2A+1
      NE1 = MAX (NE1,NE1A)
      NE2 = MIN (NE2,NE2A)
      NE2 = (NE2-NE1)/2+1
*
*   3.1  Identify type of radial integrals
*
      IBRE = 1
      IF ((IS(1) .EQ. IS(4)) .AND. (IS(2) .NE. IS(3))) IBRE = 2
      IF ((IS(1) .NE. IS(4)) .AND. (IS(2) .EQ. IS(3))) IBRE = 2
      IF ((IS(1) .EQ. IS(3)) .AND. (IS(2) .EQ. IS(4))) IBRE = 4
      RETURN
*
      END
