************************************************************************
*                                                                      *
      SUBROUTINE SETJ (IS,JS,KS,NS,KJ23)
*                                                                      *
*   Sets the tables required by  the recoupling  coefficient package   *
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
*   NJGRAF. This routine loads  the COMMON block /COUPLE/ with para-   *
*   meters for the first call  of NJGRAF involving direct integrals.   *
*   Subsequent exchange calls  of NJGRAF must be preceeded by a call   *
*   of MODJ23 to restore these arrays to their correct initial state.  *
*                                                                      *
*                                           Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      PARAMETER (MANGM = 60, MTRIAD = 12)
*
      LOGICAL FREE
*
      DIMENSION IS(2,2),JS(2,2),KS(2,2)
*
      COMMON/L1/JBQ1(3,NNNW),JBQ2(3,NNNW),JTQ1(3),JTQ2(3)
     :      /L2/J2S(MTRIAD,3),J3S(MTRIAD,3)
     :      /M0/JJC1(NNNW),JJC2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
*
*   1.0  Set J1 array
*
      II = 0
      DO 1 IJ = 1,NS
         I = JLIST(IJ)
         II = II+1
         J1(II) = JBQ1(3,I)
    1 CONTINUE
      IF (NS .EQ. 1) GOTO 4
      NS1 = NS-1
      DO 2 I = 1,NS1
         II = II+1
         J1(II) = JJC1(I)
    2 CONTINUE
      DO 3 I = 1,NS1
         II = II+1
         J1(II) = JJC2(I)
    3 CONTINUE
    4 CONTINUE
      DO 5 I = 1,2
         II = II+1
         IJ = IS(I,1)
         J1(II) = JJQ1(3,IJ)
         IF ((I .EQ. 1) .AND. (IS(1,1) .EQ. IS(2,1))) J1(II) = JTQ1(3)
         J1(II+4) = KS(I,1)
    5 CONTINUE
      DO 6 I = 1,2
         II = II+1
         IJ = IS(I,2)
         J1(II) = JJQ2(3,IJ)
         IF ((I .EQ. 1) .AND. (IS(1,2) .EQ. IS(2,2))) J1(II) = JTQ2(3)
         J1(II+4) = KS(I,2)
    6 CONTINUE
*
*   2.0  Set J2, J3 arrays if not already available
*
      NS2 = MAX (4,NS+2)
      IF (KJ23 .GT. 0) GOTO 14
*
      DO 7 I = 4,NS2
         J2(I,1) = NS+I-4
         J2(I,2) = I-2
         J2(I,3) = NS+I-3
         J3(I,1) = J2(I,1)+NS-1
         J3(I,2) = I-2
         J3(I,3) = J2(I,3)+NS-1
    7 CONTINUE
      J2(4,1) = 1
      J3(4,1) = 1
*
*   At this stage, the entries in rows corresponding to active
*   shells are set incorrectly.
*
*   3.0  Set rows 1 through 3
*
      NS3 = 3*NS
      J2(1,1) = NS3+5
      J2(1,2) = NS3+7
      J2(1,3) = NS3+3
      J2(2,1) = JS(1,1)
      J2(2,2) = NS3+3
      J2(2,3) = NS3-1
      J2(3,1) = JS(2,1)
      J2(3,2) = NS3+4
      J2(3,3) = NS3
*
      J3(1,1) = NS3+7
      J3(1,2) = NS3+4
      J3(1,3) = NS3+6
      J3(2,1) = JS(1,2)
      J3(2,2) = NS3+5
      J3(2,3) = NS3+1
      J3(3,1) = JS(2,2)
      J3(3,2) = NS3+6
      J3(3,3) = NS3+2
*
*   4.0  Set remaining resultants
*
      IJ1 = JS(1,1)
      IJ2 = JS(2,1)
      IF (IJ2 .GT. 1) J2(IJ2+2,2) = J2(3,3)
      IF (IJ2 .EQ. 1) J2(4,1) = J2(3,3)
      IF (IJ1 .NE. IJ2) GOTO 8
      J2(3,1) = J2(2,3)
      GOTO 9
*
    8 IF (IJ1 .GT. 1) J2(IJ1+2,2) = J2(2,3)
      IF (IJ1 .EQ. 1) J2(4,1) = J2(2,3)
*
    9 IJ1 = JS(1,2)
      IJ2 = JS(2,2)
      IF (IJ2 .GT. 1) J3(IJ2+2,2) = J3(3,3)
      IF (IJ2 .EQ. 1) J3(4,1) = J3(3,3)
      IF (IJ1 .NE. IJ2) GOTO 10
      J3(3,1) = J3(2,3)
      GOTO 11
*
   10 IF (IJ1 .GT. 1) J3(IJ1+2,2) = J3(2,3)
      IF (IJ1 .EQ. 1) J3(4,1) = J3(2,3)
*
*   All arrays now set. Put up flag KJ23.
*
   11 KJ23 = 1
      MJA = NS3+7
      NJA = NS+3
*
*   5.0  Save J2, J3 and return
*
      DO 13 J = 1,3
         DO 12 I = 1,NS2
            J2S(I,J) = J2(I,J)
            J3S(I,J) = J3(I,J)
   12    CONTINUE
   13 CONTINUE
      RETURN
*
*   6.0  Reset J2, J3 from buffers if KJ23 has been set
*
   14 DO 16 J = 1,3
         DO 15 I = 1,NS2
            J2(I,J) = J2S(I,J)
            J3(I,J) = J3S(I,J)
   15    CONTINUE
   16 CONTINUE
      RETURN
*
      END
