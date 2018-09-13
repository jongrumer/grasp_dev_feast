************************************************************************
*                                                                      *
      SUBROUTINE RKCO_GG (JA,JB,CORD,INCOR,ICOLBREI)
*                                                                      *
*   Configurations JA, JB. Analyse the tables of quantum numbers set   *
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
*   in the COMMON  blocks M0 , M1, M2, M3  to determine all possible   *
*   sets of interacting  orbitals which give a non-vanishing Coulomb   *
*   matrix element,  and  initiates the calculation of coefficients.   *
*   The following conventions are in force: (1) labels 1, 2 refer to   *
*   left, right sides of matrix element respectively;   (2) pointers   *
*   JA1, JB1, JA2, JB2 point to the JLIST array of active  orbitals;   *
*   IA1, IB1, IA2, IB2 point to the complete list of orbitals.         *
*                                                                      *
*   Call(s) to: [LIB92]: COR, CORD, ISPAR, ITJPO, SETQNA, VIJOUT.      *
*                                                                      *
*                                           Last update: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /DUMX/JLIS(NNNW),JC1S(NNNW),JC2S(NNNW)
     :      /M0/JJC1(NNNW),JJC2(NNNW)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
CGG      EXTERNAL COR, CORD
      EXTERNAL CORD
*
*   The Hamiltonian is an even scalar operator
*
      IF (ITJPO (JA) .NE. ITJPO (JB)) RETURN
      IF (ISPAR (JA) .NE. ISPAR (JB)) RETURN
      IF (ICHKQ2(JA,JB) .EQ. 0) RETURN
*
      CALL SETQNA (JA,JB)
      IF (IBUG2 .EQ. 1) CALL VIJOUT (JA,JB)
*
*   1.0 Analyse peel shell interactions
*
*   1.1 Analyse electron distribution in peel. (The full procedure is
*       needed only if the number of peel orbitals NPEEL .GE. 2)
*
      IF (NW .LT. 1) THEN
         PRINT *, 'RKCO: No subshells.'
         STOP
      ENDIF
      IF (NPEEL .EQ. 0) GOTO 48
      IF (NPEEL .EQ. 1) GOTO 43
*
*   Find differences in occupations, NDQ, for each peel orbital in
*   turn and use to set up labels of active orbitals maintaining the
*   convention JA1 .LE. JB1, JA2 .LE. JB2.
*
      IDQ = 0
      JA1 = 0
      JB1 = 0
      JA2 = 0
      JB2 = 0
CGG
      IDQG = 0
CGG
      DO 10 JW = 1,NPEEL
         J = JLIST(JW)
         NDQ = NQ1(J) - NQ2(J)
         IF (IABS (NDQ) .GT. 2) RETURN
         IF (NDQ .LT. 0) GOTO 5
         IF (NDQ-1) 10,1,4
    1    IF (JA1 .GT. 0) GOTO 2
         JA1 = JW
         GOTO 3
    2    JB1 = JW
    3    IDQ = IDQ+1
CGG
         IDQG=1+IDQG
CGG
         GOTO 10
    4    JA1 = JW
         IDQ = IDQ+2
CGG
         IDQG=20+IDQG
CGG
         GOTO 10
    5    IF (NDQ+1) 9,6,10
    6    IF (JA2 .GT. 0) GOTO 7
         JA2 = JW
         GOTO 8
    7    JB2 = JW
    8    IDQ = IDQ+1
CGG
         IDQG=1+IDQG
CGG
         GOTO 10
    9    JA2 = JW
         IDQ = IDQ+2
CGG
         IDQG=20+IDQG
CGG
   10 CONTINUE
*
*   1.2 Calculate coefficients for all possible sets of active shells.
*
*   There are 4 cases, depending on the value of IDQ, the sum of the
*   absolute differences NDQ:
*
*   1.2.1 IDQ .GT. 4: matrix element null
*
      IF (IDQ .GT. 4) RETURN
      IF (IDQ .EQ. 4) GOTO 12
      IF (IDQ .NE. 2) GOTO 11
      KLAST = 1
      GOTO 16
   11 IF (IDQ .NE. 0) GOTO 54
      IF (JA .EQ. JB) GOTO 43
      KLAST = NPEEL
      GOTO 16
*
*   1.2.2 IDQ .EQ. 4: matrix element uniquely defined
*
   12 IF (JB1 .NE. 0) GOTO 13
      JB1 = JA1
   13 IF (JB2 .NE. 0) GOTO 14
      JB2 = JA2
   14 IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA2,JB2
CGG
CGG      CALL COR (JA,JB,JA1,JB1,JA2,JB2)
      IF(IDQG.NE.40)GO TO 1003
C
C     TARP KONFIGURACIJU
C                        KAI      N'1=N1+1
C                                 N'2=N2-1
C
      CALL EL2(JA,JB,JA1,JA2,ICOLBREI)
      RETURN
 1003 IF(IDQG.NE.22)GO TO 1009
      CALL EL4(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
      RETURN
C
C     TARP KONFIGURACIJU
C                        KAI      N'1=N1+-1
C                                 N'2=N2+-1
C                        KAI      N'3=N3+-1
C                                 N'4=N4+-1
C
 1009 CALL EL5(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
      RETURN
CGG
*
*   1.2.3 IDQ .EQ. 2: One orbital fixed each side include all
*                     possible spectators.
*
*   Also IDQ .EQ. 0 for a matrix element off-diagonal in coupling
*   only. Must sum over all pairs of orbitals excluding core-core
*   terms
*
   16 DO 42 KWA = 1,KLAST
         IF (IDQ .EQ. 2) GOTO 17
         JA1 = KWA
         JA2 = KWA
   17    JT1 = JA1
         JT2 = JA2
         IT1 = JLIST(JA1)
         IT2 = JLIST(JA2)
         DO 25 KW = KWA,NPEEL
            K1 = JLIST(KW)
            IF (NQ1(K1)*NQ2(K1) .EQ. 0) GOTO 25
            JB1 = KW
            JB2 = KW
            JA1 = JT1
            JA2 = JT2
*
*   Interchange JA1 and JB1 and/or JA2 and JB2 if necessary
*
            IF (JA1-JB1) 20,19,18
   18       JT3 = JB1
            JB1 = JA1
            JA1 = JT3
            GOTO 20
   19       IB1 = JLIST(JB1)
            IF (NQ1(IB1) .LE. 1) GOTO 25
   20       IF (JA2-JB2) 23,22,21
   21       JT3 = JB2
            JB2 = JA2
            JA2 = JT3
            GOTO 23
   22       IB2 = JLIST(JB2)
            IF (NQ2(IB2) .LE. 1) GOTO 25
   23       IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA2,JB2
CGG
CGG            CALL COR (JA,JB,JA1,JB1,JA2,JB2)
             IF(IDQ.NE.0)GO TO 1002
C
C     TARP TU PACIU KONFIGURACIJU
C
             CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
             GO TO 25
C
C     TARP KONFIGURACIJU
C                        KAI      N'1=N1+1
C                                 N'2=N2-1
C                                 N'3=N3
C
 1002        CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
CGG
   25    CONTINUE
         IF ((IDQ .EQ. 0) .AND. (NCORE .EQ. 0)) GOTO 42
         IF ((NCORE .EQ. 0) .OR. (NAK(IT1) .NE. NAK(IT2))) RETURN
*
*   This section calculates the terms arising from active electrons
*   which are in closed shells
*
         NPEELM = NPEEL-1
         DO 26 I = 1,NPEEL
            JLIS(I) = JLIST(I)
   26    CONTINUE
         DO 27 I = 1,NPEELM
            JC1S(I) = JJC1(I)
            JC2S(I) = JJC2(I)
   27    CONTINUE
         DO 41 KW = 1,NCORE
            IJW = KLIST(KW)
            DO 28 I = 1,NPEEL
               IJ = JLIST(I)
               IF (IJW .LT. IJ) GOTO 29
   28       CONTINUE
            I = NPEEL+1
            GOTO 31
   29       IM = NPEEL-I+1
            DO 30 II = 1,IM
               JLIST(NPEEL+2-II) = JLIST(NPEEL+1-II)
               IF (NPEEL .EQ. II) GOTO 31
               JJC1(NPEEL+1-II) = JJC1(NPEEL-II)
               JJC2(NPEEL+1-II) = JJC2(NPEEL-II)
   30       CONTINUE
   31       CONTINUE
            IF (I .LT. 3) GOTO 32
            JJC1(I-1) = JJC1(I-2)
            JJC2(I-1) = JJC2(I-2)
            GOTO 33
   32       I1 = JLIST(1)
            JJC1(1) = JJQ1(3,I1)
            JJC2(1) = JJQ2(3,I1)
   33       JLIST(I) = IJW
            JA1 = JT1
            IF (JT1 .GE. I) JA1 = JA1+1
            JB1 = I
            JA2 = JT2
            IF (JT2 .GE. I) JA2 = JA2+1
            JB2 = I
            IF (JA1-JB1) 35,35,34
   34       JT3 = JB1
            JB1 = JA1
            JA1 = JT3
   35       CONTINUE
            IF (JA2-JB2) 37,37,36
   36       JT3 = JB2
            JB2 = JA2
            JA2 = JT3
   37       CONTINUE
            NPEEL = NPEEL+1
            IF (IBUG2 .NE. 0) THEN
               NPEELM = NPEEL-1
               WRITE (99,302) JA1,JB1,JA2,JB2,KW,KLIST(KW)
               WRITE (99,303) (JLIST(I),I = 1,NPEEL)
               WRITE (99,304) (JJC1(I),I = 1,NPEELM)
               WRITE (99,305) (JJC2(I),I = 1,NPEELM)
            ENDIF
CGG
CGG            CALL COR (JA,JB,JA1,JB1,JA2,JB2)
      IF(IDQ.NE.0)GO TO 1001
C
C     TARP TU PACIU KONFIGURACIJU
C
      CALL EL1(JA,JB,JA1,JB1,1,ICOLBREI)
      GO TO 1146
 1001 IF(IDQG.NE.40)GO TO 1005
 1007 WRITE(99,995)
  995 FORMAT('   rie zymes 38 atv N=N-N  !!!!!!')
      RETURN
 1005 IF(IDQG.NE.2)GO TO 1007
C
C     TARP KONFIGURACIJU
C                        KAI      N'1=N1+1
C                                 N'2=N2-1
C                                 N'3=N3
C
       CALL EL3(JA,JB,JA1,JB1,JA2,JB2,ICOLBREI)
 1146       NPEEL = NPEEL-1
CGG
CGG            NPEEL = NPEEL-1
            NPEELM = NPEEL-1
            DO 39 I = 1,NPEEL
               JLIST(I) = JLIS(I)
   39       CONTINUE
            DO 40 I = 1,NPEELM
               JJC1(I)  = JC1S(I)
               JJC2(I)  = JC2S(I)
   40       CONTINUE
   41    CONTINUE
   42 CONTINUE
      RETURN
*
*   1.2.4 IDQ .EQ. 0 - diagonal case. Include all pairs with
*         JA1 = JA2, JB1 = JB2.
*
   43 DO 47 KW1 = 1,NPEEL
         K1 = JLIST(KW1)
         JB1 = KW1
         JB2 = KW1
         DO 46 KW2 = 1,KW1
            JA1 = KW2
            IF (JA1 .NE. JB1) GOTO 44
            IF (NQ1(K1) .LE. 1) GOTO 46
   44       JA2 = JA1
            IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA2,JB2
CGG
CGG            CALL COR (JA,JB,JA1,JB1,JA2,JB2)
      IF(JA.NE.JB)GO TO 1000
C
C     TARP TU PACIU BUSENU
C
      CALL EL1(JA,JB,JA1,JB1,0,ICOLBREI)
      GO TO 46
 1010 WRITE(99,996)
  996 FORMAT('   rie zymes 45 atv N=N-N  !!!!!!')
      WRITE(99,*)"JA,JB,JA1,JB1",JA,JB,JA1,JB1
      WRITE(99,*)"IDQG IDQ",IDQG,IDQ
      RETURN
 1000 IF(IDQG.NE.2)GO TO 1010
C
C     TARP KONFIGURACIJU
C                        KAI      N'1=N1+1
C                                 N'2=N2-1
C                                 N'3=N3
C
       CALL EL3(JA1,JB1,JA2,JB2,ICOLBREI)
   46    CONTINUE
   47 CONTINUE
   48 IF (INCOR .LT. 1) RETURN
      IF (NCORE .EQ. 0) RETURN
*
*   2.0 The diagonal case. deal with contributions from core orbitals
*       if INCOR .EQ. 1.
*
CGG
CGG  Cia orginalioje programoje yra klaida
CGG
      IF(JA.NE.JB) RETURN
CGG
      DO 53 KW1 = 1,NCORE
         JB1 = KW1
         JB2 = KW1
*
*   2.1 Calculate contribution from core/core terms
*
         IPCA = 2
         DO 50 KW2 = 1,KW1
            JA1 = KW2
            JA2 = KW2
            IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA1,JB1
            CALL CORD (JA,JB,JA1,IPCA,JB1)
   50    CONTINUE
*
*   2.2 Calculate contribution from peel/core terms
*
         IF (NPEEL .EQ. 0) GOTO 53
         IPCA = 1
         DO 52 KW2 = 1,NPEEL
            JA1 = KW2
            JA2 = KW2
            IF (IBUG2 .NE. 0) WRITE (99,301) JA1,JB1,JA1,JB1
            CALL CORD (JA,JB,JA1,IPCA,JB1)
   52    CONTINUE
   53 CONTINUE
      RETURN
*
*   3.0 Diagnostic print - NW .LT. 1
*
   54 WRITE (*,300)
      STOP
*
  300 FORMAT ('RKCO: Error.')
  301 FORMAT ('From RKCO:'
     :         /10X,' JA1 = ',I3,4X,' JB1 = ',I3,4X,' JA2 = ',I3,
     :              ' JB2 = ',I3)
  302 FORMAT ('From RKCO:'
     :   /10X,' JA1 = ',I3,4X,' JB1 = ',I3,4X,' JA2 = ',I3,
     :        ' JB2 = ',I3,' K2  = ',I3,   ' KW  = ',I3)
  303 FORMAT (1X,'JLIST : ',25I4)
  304 FORMAT (1X,'JJC1  : ',25I4)
  305 FORMAT (1X,'JJC2  : ',25I4)
*
      END
