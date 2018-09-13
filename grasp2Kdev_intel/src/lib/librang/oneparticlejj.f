************************************************************************
*                                                                      *
      SUBROUTINE ONEPARTICLEJJ(KA,IOPAR,JA,JB,IA1,IA2,VSHELL)
*                                                                      *
*   The  main  program for evaluating the reduced matrix elements of   *
*   a one particle operator for configurations in jj-coupling.         *
*                                                                      *
*   Call(s) to: [LIB92]: CFP, FIXJ, GENSUM, ICHOP, IROW1, ISPAR,       *
*                        ITJPO, ITRIG, SETQNA, VIJOUT.                 *
*               [NJGRAF]: NJGRAF.                                      *
*                                                                      *
*                                           Last update: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      PARAMETER (MANGM = 60, MTRIAD = 12)
      PARAMETER (M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3))
      PARAMETER (M6J = 20, MSUM = 10)
*
      LOGICAL FREE,SUMVAR,FAIL
*
      DIMENSION VSHELL(NNNW)
      DIMENSION IS(2),KS(2)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
*
      COMMON/ARGU/J6C,J7C,J8C,J9C,JWC,J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :       J9(MANGMP),KW(6,M6J),JDEL,LDEL(M6J,2),SUMVAR(MANGM),MP
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DUMX/JLIS(NNNW),JC1S(NNNW),JC2S(NNNW)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /L1/JBQ1(3,NNNW),JBQ2(3,NNNW),JTQ1(3),JTQ2(3)
     :      /M0/JJC1(NNNW),JJC2(NNNW)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
      COMMON/SUMARG/J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :              JWORD(6,M6J),NLSUM,NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :              K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :              JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
      IA1 = 0
      KK = KA+KA+1
      IF (ITRIG (ITJPO (JA),ITJPO (JB),KK) .EQ. 0) RETURN
      IF ((IOPAR .NE. 0) .AND. (ISPAR (JA)*ISPAR (JB)*IOPAR .NE. 1))
     :   RETURN
      IF (ICHKQ1(JA,JB).EQ. 0) RETURN
*
      CALL SETQNA (JA,JB)
      IF (IBUG4 .NE. 0) CALL VIJOUT (JA,JB)
*
      DO 1 IJ = 1,NW
         VSHELL(IJ) = ZERO
    1 CONTINUE
*
*   Analyse peel shell interactions
*
      IDQ = 0
      JA1 = 0
      JA2 = 0
*
      IF (NPEEL .NE. 0) THEN
*
        DO 2 JW = 1,NPEEL
          IJ = JLIST(JW)
          NDQ = NQ1(IJ)-NQ2(IJ)
          IF (ABS (NDQ) .GT. 1) GOTO 39
          IF (NDQ .GT. 0) THEN
            JA1 = JW
            IDQ = IDQ+1
          ELSEIF (NDQ .LT. 0) THEN
            JA2 = JW
            IDQ = IDQ+1
          ENDIF
    2   CONTINUE
*
        IF (IDQ .GT. 2) GOTO 39
*
*   Evaluate the array VSHELL
*
*   Then there are two possibilities IDQ = 0 or IDQ = 2
*   if IDQ = 0, then loop over all shells by index ISH
*   if IDQ = 2, then one orbital fixed on each side
*
        NS = NPEEL
      ENDIF
*
      IF (IDQ .EQ. 2) GOTO 19
*
*   Loop over shells when IDQ = 0
*
      ISH = 0
      IF (NPEEL .EQ. 0) GOTO 9
      DO 7 I = 1,NPEEL
    7    JLIS(I) = JLIST(I)
      IF (NPEEL .EQ. 1) GOTO 9
      NPEELM = NPEEL-1
      DO 8 I = 1,NPEELM
         JC1S(I) = JJC1(I)
    8    JC2S(I) = JJC2(I)
*
*   If ISH .GT. NW, then loop is over and return
*
    9 ISH = ISH+1
      IF (ISH .GT. NW) RETURN
      IF (ICHOP (ISH,JA) .EQ. -1) GOTO 9
      IF (IBUG6 .NE. 0) WRITE (99,308) ISH
      IF (ICHOP (ISH,JA) .EQ. 0) GOTO 16
*
*   Case one --- the ISH-th shell is in the core or in the peel and
*   closed for both sides
*
      I = 1
      IF (NPEEL.EQ.0) GOTO 15
      DO 10 I = 1,NPEEL
        IJ = JLIST(I)
        IF (ISH .LT. IJ) GOTO 11
   10 CONTINUE
      I = NPEEL+1
      GOTO 13
   11 IM = NPEEL-I+1
      DO 12 II = 1,IM
         JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
         IF (NPEEL.EQ.II) GOTO 13
         JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
         JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
   12 CONTINUE
   13 CONTINUE
      IF (I .LT. 3) GOTO 14
      JJC1(I-1) = JJC1(I-2)
      JJC2(I-1) = JJC2(I-2)
      GOTO 15
   14 I1 = JLIST(1)
      JJC1(1) = JJQ1(3,I1)
      JJC2(1) = JJQ2(3,I1)
   15 JLIST(I) = ISH
      JA1 = I
      JA2 = I
      NS = NPEEL+1
      GOTO 19
*
*   Case two --- the ISH-th shell is in the peel and open for either
*   side
*
   16 NS = NPEEL
      DO 17 JW = 1,NPEEL
        NX = ISH-JLIST(JW)
        IF (NX.EQ.0) GOTO 18
   17 CONTINUE
   18 JA1 = JW
      JA2 = JW
*
*   Main computation
*
*     JA1, JA2 are the indices of interacting shells in JLIST
*     IA1, IA2 are the indices of interacting shells in NW
*
   19 IA1 = JLIST(JA1)
      IA2 = JLIST(JA2)
      KS1 = 2*ABS (NAK(IA1))
      KS2 = 2*ABS (NAK(IA2))
*
*   Check triangular condition for the active shells
*
      IF (ITRIG (KS1,KS2,KK).EQ.1) GOTO 99
      IF (IDQ .EQ. 2) RETURN
      GOTO 100
*
*   Set tables of quantum numbers of non-interacting spectator shells
*
   99 CONTINUE
      IF (IDQ .EQ. 0) THEN
CGG
CGG   Cia orginalioje programoje yra klaida
CGG
CGG        IF(JA .EQ. JB) THEN
          CALL ONEPARTICLEJJ1(NS,KA,JA,JB,JA1,JA2,TCOEFF)
CGG          CALL ONESCALAR1(NS,JA,JB,JA1,JA2,TCOEFF)
CGG        ELSE
CGG          TCOEFF = 0.0D 00 
CGG        END IF
      ELSE IF (IDQ .EQ. 2) THEN 
*
*   IDQ = 2 Case
*
*       Permutation factor for IDQ = 2
*
        CALL ONEPARTICLEJJ2(NS,KA,JA1,JA2,TCOEFF)
CGG        CALL ONESCALAR2(JA,JB,JA1,JA2,TCOEFF)
        VSHELL(1) = TCOEFF 
        RETURN
      END IF
*
*   End of loop over parent states
*
*
*   IDQ = 0 CASE
*
      VSHELL(ISH) = TCOEFF
*
*   Loop over all shells when IDQ = 0
*
100   CONTINUE
      IF (NPEEL .EQ. 0) GOTO 9
      DO 35 I = 1,NPEEL
   35    JLIST(I) = JLIS(I)
      IF (NPEEL .EQ. 1) GOTO 9
      NPEELM = NPEEL-1
      DO 36 I = 1,NPEELM
         JJC1(I)  = JC1S(I)
   36    JJC2(I)  = JC2S(I)
      GOTO 9
*
   39 IF (IBUG6 .NE. 0) WRITE (99,300)
      RETURN
   40 IF (IBUG6 .NE. 0) WRITE (99,301)
      RETURN
*
  300 FORMAT (' One side has more than one interacting electron')
  301 FORMAT (' Spectator quantum numbers not diagonal for non-interact'
     :   ,'ing shells')
  302 FORMAT (/' J1')
  303 FORMAT (24I5)
  304 FORMAT (' J2                   J3')
  305 FORMAT (3I5,I10,2I5)
  306 FORMAT(' CFP  ',I3,I4,I7,2I4,I7,2I4,1P,D20.9)
  307 FORMAT(/' Recoupling coefficient = ',1P,D19.12)
  308 FORMAT(//' ISH = ',I3)
*
      END
