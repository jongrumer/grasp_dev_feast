************************************************************************
*                                                                      *
      SUBROUTINE FIXJ (JA1,JA2,KA,IS,KS,NS,KJ23)
*                                                                      *
*   Sets up the arrays J1, J2, J3 required by the recoupling package   *
*   NJSYM.                                                             *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      PARAMETER (MANGM = 60, MTRIAD = 12)
*
      LOGICAL FREE
*
      DIMENSION IS(2),KS(2)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/M0/JJC1(NNNW),JJC2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /L1/JBQ1(3,NNNW),JBQ2(3,NNNW),JTQ1(3),JTQ2(3)
*
*   Set up the J2 and J3 arrays
*
      NM1 = NS-1
      IF (KJ23 .EQ. 1) GOTO 5
      NS1 = NS+1
      N2 = NS+NS
      N3 = N2+NS
*
      J2(1,1) = N3+2
      J2(1,2) = N3+3
      J2(1,3) = N3+1
      J2(2,1) = JA1
      J2(2,2) = N3+1
      J2(2,3) = N3-1
*
      J3(1,1) = JA2
      J3(1,2) = N3+2
      J3(1,3) = N3
*
      IF (NS .EQ. 1) GOTO 3
*
      DO 1 JW = 1,NM1
         JJ = JW+2
         J2(JJ,1) = NS+JW-1
         J2(JJ,2) = JW+1
         J2(JJ,3) = NS+JW
*
         JK = JW+1
         J3(JK,1) = N2+JW-2
         J3(JK,2) = JW+1
         J3(JK,3) = N2+JW-1
    1 CONTINUE
*
      J2(3,1) = 1
      IF (JA1 .EQ. 1) J2(3,1) = N3-1
*
      J3(2,1) = 1
      IF (JA2 .EQ. 1) J3(2,1) = N3
*
      J2(NS1,3) = N2-1
*
      J3(NS1,1) = N3-2
      J3(NS1,2) = N3+3
      J3(NS1,3) = N2-1
*
      IF (JA1 .EQ. 1) GOTO 2
      JAF1 = JA1+1
      J2(JAF1,2) = N3-1
*
    2 IF (JA2 .EQ. 1) GOTO 4
      J3(JA2,2) = N3
*
      IF (NS .GT. 1) GOTO 4
    3 J3(2,1) = N3
      J3(2,2) = N3+3
      J3(2,3) = N3-1
*
    4 CONTINUE
*
*   Set the J1 array
*
    5 CONTINUE
      II = 0
*
      DO 6 JW = 1,NS
         IJ = JLIST(JW)
         II = II+1
         J1(II) = JBQ2(3,IJ)
    6 CONTINUE
*
      IF (NS .EQ. 1) GOTO 9
*
      DO 7 JW = 1,NM1
         II = II+1
         J1(II) = JJC1(JW)
    7 CONTINUE
*
      DO 8 JW = 1,NM1
         II = II+1
         J1(II) = JJC2(JW)
    8 CONTINUE
*
    9 CONTINUE
      II = II+1
      IJ = IS(1)
      J1(II) = JJQ1(3,IJ)
      J1(II+2) = KS(1)
      II = II+1
      IJ = IS(2)
      J1(II) = JJQ2(3,IJ)
      J1(II+2) = KS(2)
*
      II = II+3
      J1(II) = KA+KA+1
      MJA = II
      NJA = NS+2
*
      RETURN
*
      END
