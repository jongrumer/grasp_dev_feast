************************************************************************
*                                                                      *
      SUBROUTINE COR (JA,JB,JA1,JB1,JA2,JB2)
*                                                                      *
*   Computes  the  MCP  coefficients.  Equation numbers are those of   *
*   Computer Phys Commun 5 (1973) 263                                  *
*                                                                      *
*   Call(s) to: [LIB92]: CRE, ITRIG, LTAB, MODJ23, MUMDAD, OCON,       *
*                        SETJ, SKRC, SPEAK, KNJ.                       *
*               [NJGRAF]: NJGRAF, GENSUM.                              *
*                                                                      *
*                                           Last update: 02 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      PARAMETER (
     :   MANGM=60,M3MNGM=3*MANGM,MANGMP=2*(MANGM/3),
     :   MTRIAD=12,
     :   M6J=20,MSUM=10)
*
      PARAMETER (EPS = 1.0D-10)
*
      LOGICAL FREE,DSUMVR,ESUMVR,FAILD,FAILE
*
      PARAMETER (IDIM = 11)
      DIMENSION COND(IDIM),CONE(IDIM)
*
      DIMENSION KAPS(4),KS(4),NQS(4),ILS(4),LLS(4),IROWS(4)
      DIMENSION IS(4),JS(4)
*
      DIMENSION JD6(M3MNGM),JD7(M3MNGM),JD8(M3MNGM),
     : JD9(MANGMP),KDW(6,M6J),LDDEL(M6J,2),DSUMVR(MANGM)
      DIMENSION JD6P(MANGMP),JD7P(MANGMP),JD8P(MANGMP),JD9P(MANGMP),
     : JDWORD(6,M6J),
     : NDBJ(MSUM),NDB6J(MSUM),KD6CP(MSUM),KD7CP(MSUM),KD8CP(MSUM),
     : KD9CP(MSUM),JDSUM6(MTRIAD),JDSUM4(MTRIAD,M6J),JDSUM5(MTRIAD,M6J),
     : INVD6J(M6J)
      DIMENSION JE6(M3MNGM),JE7(M3MNGM),JE8(M3MNGM),
     : JE9(MANGMP),KEW(6,M6J),LEDEL(M6J,2),ESUMVR(MANGM)
      DIMENSION JE6P(MANGMP),JE7P(MANGMP),JE8P(MANGMP),JE9P(MANGMP),
     : JEWORD(6,M6J),
     : NEBJ(MSUM),NEB6J(MSUM),KE6CP(MSUM),KE7CP(MSUM),KE8CP(MSUM),
     : KE9CP(MSUM),JESUM6(MTRIAD),JESUM4(MTRIAD,M6J),JESUM5(MTRIAD,M6J),
     : INVE6J(M6J)
*
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /COUPLE/MJA,NJA,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /L1/JBQ1(3,NNNW),JBQ2(3,NNNW),JTQ1(3),JTQ2(3)
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
*   1.0  Initialize pointers and flags and set any
*        tables required.
*
*        In this segment, the array IS points to the
*        full list of orbitals, the array JS to the
*        array JLIST of peel orbital pointers.
*
*   1.1  Initialization
*
      JS(1) = JA1
      JS(2) = JB1
      JS(3) = JA2
      JS(4) = JB2
      DO 1 I = 1,4
        IS(I) = JLIST(JS(I))
        KAPS(I) = 2*NAK(IS(I))
        KS(I) = ABS (KAPS(I))
    1 CONTINUE
      IA1 = IS(1)
      IB1 = IS(2)
      IA2 = IS(3)
      IB2 = IS(4)
*
      KJ23 = 0
      ISNJ = 0
      FAILD = .FALSE.
      FAILE = .FALSE.
*
*   Initialize arrays
*
      DO 3 J = 1,NW
         DO 2 K = 1,3
            JBQ1(K,J) = 0
            JBQ2(K,J) = 0
    2    CONTINUE
    3 CONTINUE
*
      NBRJ = 3*NPEEL+7
      DO 4 I = 1,NBRJ
         FREE(I) = .FALSE.
    4 CONTINUE
*
*   2.0 Set tables of quantum numbers of spectator
*       shells.
*
      DO 8 JJ = 1,NPEEL
        J = JLIST(JJ)
        IF ((J .NE. IA1) .AND. (J .NE. IB1)) THEN
          DO 5 K = 1,3
            JBQ1(K,J) = JJQ1(K,J)
    5     CONTINUE
        ENDIF
        IF ((J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 6 K = 1,3
            JBQ2(K,J) = JJQ2(K,J)
    6     CONTINUE
        ENDIF
*
*   2.1 Examine quantum numbers of spectator
*       shells for orthogonality and
*       exit if found.
*
        IF ((J .NE. IA1) .AND. (J .NE. IB1) .AND.
     :      (J .NE. IA2) .AND. (J .NE. IB2)) THEN
          DO 7 K = 1,3
            IF (JBQ1(K,J) .NE. JBQ2(K,J)) THEN
               IF (IBUG2 .NE. 0) WRITE (99,300)
               GOTO 41
            ENDIF
    7     CONTINUE
        ENDIF
    8 CONTINUE
*
*   3.1  Set range of the parameter k for Coulomb
*        integrals.
*        Terminate run if buffer store dimension
*        IDIM is too small.
*
      CALL SKRC (IS,KAPS,KS,KD1,KD2,KE1,KE2)
      IF ((KD2 .EQ. 0) .AND. (KE2 .EQ. 0)) GOTO 41
      IF ((KD2. GT. IDIM) .OR. (KE2 .GT. IDIM)) THEN
         KK = MAX (KE2,KD2)
         WRITE (*,301) KK
         STOP
      ENDIF
*
      IF (KD2 .NE. 0) THEN
         DO 9 K = 1,KD2
            COND(K) = ZERO
    9    CONTINUE
      ENDIF
      IF (KE2 .NE. 0) THEN
         DO 10 K = 1,KE2
            CONE(K) = ZERO
   10    CONTINUE
      ENDIF
*
      NQS(1) = NQ1(IA1)
      NQS(2) = NQ1(IB1)
      NQS(3) = NQ2(IA2)
      NQS(4) = NQ2(IB2)
*
*   3.3  Set parameters of summation over parent
*        (barred) terms in Eq. (5). The array IROWS
*        is formed to point at the list of allowed
*        parents of active shells in the array
*        NTAB.
*
      CALL LTAB (IS,NQS,KS,IROWS)
*
      DO 17 I = 1,4
         II = IROWS(I)
         LLS(I) = ITAB(II)
         ILS(I) = JTAB(II)
   17 CONTINUE
*
*   4.0  Sum contributions over all parent terms
*        permitted by angular momentum and seniority
*        selection rules.
*
      LLS1 = LLS(1)
      IF (LLS1 .NE. 1) FREE (JA1) = .TRUE.
      LLS2 = LLS(2)
      LLS3 = LLS(3)
      LLS4 = LLS(4)
*
      LS2 = ILS(2)
      DO 38 LB1 = 1,LLS2
        LS2 = LS2+3
        IT12 = NTAB(LS2)
        IT2 = KS(2)
        IT3 = JJQ1(3,IB1)
        IF (ITRIG (IT12,IT2,IT3) .EQ. 0) GOTO 38
        IF (ABS (NTAB(LS2-2)-JJQ1(1,IB1)) .NE. 1) GOTO 38
*
        LS1 = ILS(1)
        DO 37 LA1 = 1,LLS1
          LS1 = LS1+3
          IT11 = NTAB(LS1)
          IT2 = KS(1)
*
          IF (IA1 .EQ. IB1) THEN
*
*   Treat IA1 .EQ. IB1 as special case.
*
            IT3 = IT12
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 37
            IF (ABS (NTAB(LS1-2)-NTAB(LS2-2)) .NE. 1) GOTO 37
            IF (LLS2 .NE. 1) FREE(NBRJ-8) = .TRUE.
          ELSE
            IT3 = JJQ1(3,IA1)
            IF (ITRIG (IT11,IT2,IT3) .EQ. 0) GOTO 37
            IF (ABS (NTAB(LS1-2)-JJQ1(1,IA1)) .NE. 1) GOTO 37
            IF (LLS2 .NE. 1) FREE(JB1) = .TRUE.
          ENDIF
*
          LS4 = ILS(4)
          DO 36 LB2 = 1,LLS4
            LS4 = LS4+3
            IT14 = NTAB(LS4)
            IT2 = KS(4)
            IT3 = JJQ2(3,IB2)
            IF (ITRIG (IT14,IT2,IT3) .EQ. 0) GOTO 36
            IF (ABS (NTAB(LS4-2)-JJQ2(1,IB2)) .NE. 1) GOTO 36
*
            LS3 = ILS(3)
            DO 35 LA2 = 1,LLS3
              LS3 = LS3+3
              IT13 = NTAB(LS3)
              IT2 = KS(3)
*
              IF (IA2 .EQ. IB2) THEN
*
*   Treat IA2 .EQ. IB2 as special case.
*
                IT3 = IT14
                IF (LLS4 .NE. 1) FREE(NBRJ-6) = .TRUE.
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 35
                IF (ABS (NTAB(LS3-2)-NTAB(LS4-2)) .NE. 1) GOTO 35
              ELSE
                IT3 = JJQ2(3,IA2)
                IF (ITRIG (IT13,IT2,IT3) .EQ. 0) GOTO 35
                IF (ABS (NTAB(LS3-2)-JJQ2(1,IA2)) .NE. 1) GOTO 35
              ENDIF
*
*   At this point the current parent has been completely defined,
*   and its quantum numbers can now be set. The JTQ arrays must
*   be set if IA1 .EQ. IB1 or IA2 .EQ. IB2. The matrix element
*   should be diagonal in barred quantum numbers.
*
              DO 27 K = 1,3
                JBQ1(K,IA1) = NTAB(LS1+K-3)
                JBQ2(K,IA2) = NTAB(LS3+K-3)
                JTQ1(K) = 0
                IF (IB1 .EQ. IA1) THEN
                  JTQ1(K) = NTAB(LS2+K-3)
                ELSE
                  JBQ1(K,IB1) = NTAB(LS2+K-3)
                ENDIF
                JTQ2(K) = 0
                IF (IB2 .EQ. IA2) THEN
                  JTQ2(K) = NTAB(LS4+K-3)
                ELSE
                  JBQ2(K,IB2) = NTAB(LS4+K-3)
                ENDIF
                DO 26 KK = 1,4
                  IF (JBQ1(K,IS(KK)) .NE. JBQ2(K,IS(KK))) GOTO 35
   26           CONTINUE
   27         CONTINUE
*
*   4.1  Evaluate product of 4 CFPs, Eq. (5).
*
              CALL MUMDAD (IS,KAPS,PROD)
              IF (ABS (PROD) .LT. EPS) GOTO 35
*
*   4.2  Set arrays for defining the recoupling
*        coefficient.
*
              CALL SETJ (IS,JS,KS,NPEEL,KJ23)
*
              IF (ISNJ .EQ. 0) THEN
*
********************* N J G R A F   V e r s i o n **********************
*
*   Set up the arrays and variables for the direct case.
*   J1(NBRJ) ( = J1(MJA) ) is set to (2*KD1+1) so that NJGRAF is
*   called correctly.
*
                IF (KD2 .NE. 0) THEN
                  IF (KD2 .GT. 1) FREE(NBRJ) = .TRUE.
                  J1(NBRJ) = KD1+KD1+1
                  CALL NJGRAF (RECUP,-1,FAILD)
                  ISNJ = 1
                  IF (.NOT. FAILD) THEN
                  CALL KNJ (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,
     :                     KDW,JDDEL,LDDEL,DSUMVR,MDP,
     :                     JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,
     :                     NDB6J,KD6CP,KD7CP,KD8CP,KD9CP,
     :                     JDSUM4,JDSUM5,JDSUM6,INVD6J)
                  ENDIF
                ENDIF
*
*
*   Set up the arrays and variables for the exchange case.
*   J1(NBRJ) ( = J1(MJA) ) is set to (2*KE1+1) so that NJGRAF is
*   called correctly.
*
                IF (KE2 .NE. 0) THEN
                  CALL MODJ23
                  FREE(NBRJ) = .FALSE.
                  IF (KE2 .GT. 1) FREE(NBRJ) = .TRUE.
                  J1(NBRJ) = KE1+KE1+1
                  CALL NJGRAF (RECUP,-1,FAILE)
                  ISNJ = 2
                  IF (.NOT. FAILE) THEN
                  CALL KNJ (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,JE9,
     :                     KEW,JEDEL,LEDEL,ESUMVR,MEP,
     :                     JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,
     :                     NEB6J,KE6CP,KE7CP,KE8CP,KE9CP,
     :                     JESUM4,JESUM5,JESUM6,INVE6J)
                  ENDIF
                ENDIF
              ENDIF
*
*   4.3  Calculate AD, Eq. (6),
*        without the phase factor.
*
              IF ((KD2 .NE. 0) .AND. (.NOT. FAILD)) THEN
                KK = KD1-2
                DO 30 K = 1,KD2
                  KK = KK+2
                  J1(MJA) = KK+KK+1
                  CALL GENSUM (JD6C,JD7C,JD8C,JD9C,JDWC,JD6,JD7,JD8,JD9,
     +            KDW,JDDEL,LDDEL,DSUMVR,MDP,
     +            JD6P,JD7P,JD8P,JD9P,JDWORD,NDLSUM,NDBJ,NDB6J,
     +            KD6CP,KD7CP,KD8CP,KD9CP,JDSUM4,JDSUM5,JDSUM6,INVD6J,
     +            X)
                  IF (ABS (X) .GE. EPS) COND(K)=COND(K)+X*PROD
   30           CONTINUE
              ENDIF
*
*   4.4  Calculate AE, Eq. (6),
*        without the phase factor.
*
              IF ((KE2 .NE. 0) .AND. (.NOT. FAILE)) THEN
                KK = KE1-2
                DO 34 K = 1,KE2
                  KK = KK+2
                  J1(MJA) = KK+KK+1
                  CALL GENSUM (JE6C,JE7C,JE8C,JE9C,JEWC,JE6,JE7,JE8,JE9,
     +            KEW,JEDEL,LEDEL,ESUMVR,MEP,
     +            JE6P,JE7P,JE8P,JE9P,JEWORD,NELSUM,NEBJ,NEB6J,
     +            KE6CP,KE7CP,KE8CP,KE9CP,JESUM4,JESUM5,JESUM6,INVE6J,
     +            Y)
                  IF (ABS(Y) .GE. EPS) CONE(K) = CONE(K)+Y*PROD
   34           CONTINUE
              ENDIF
*
              IF (FAILD .AND. FAILE) GOTO 500
*
   35       CONTINUE
   36     CONTINUE
   37   CONTINUE
   38 CONTINUE
*
*   4.5  Insert factors independent of barred
*        quantum numbers.
*        Output results
*
*        Begin with common statistical factors, Eq. (5).
*
  500 CONST = OCON (IA1,IB1,IA2,IB2)
*
      KAP1 = NAK(IA1)
      KAP2 = NAK(IB1)
      KAP3 = NAK(IA2)
      KAP4 = NAK(IB2)
*
*   4.6  Compute products of reduced matrix
*        elements, Eq. (7).
*        CRED for direct terms
*        CREE for exchange terms
*
      IF (KD2 .NE. 0) THEN
         PRODD = CONST/SQRT (DBLE (KS(1)*KS(4)))
         IF (MOD (KD1,2) .NE. 0) PRODD = -PRODD
         IF ((IA1 .EQ. IB1) .AND. (IA2 .EQ. IB2)) PRODD = PRODD*HALF
         KK = KD1-2
         DO 39 K = 1,KD2
            KK = KK+2
            CRED = CRE (KAP1,KK,KAP3)*CRE (KAP2,KK,KAP4)
            X = PRODD*COND(K)*CRED
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IA2,IB2,KK,X)
   39    CONTINUE
      ENDIF
*
      IF (KE2 .NE. 0) THEN
         PRODE = CONST/SQRT(DBLE (KS(1)*KS(3)))
         IF (MOD (KE1,2) .NE. 0) PRODE = -PRODE
         PRODE = -PRODE
         KK = KE1-2
         DO 40 K = 1,KE2
            KK = KK+2
            CREE = CRE (KAP1,KK,KAP4)*CRE (KAP2,KK,KAP3)
            X = PRODE*CONE(K)*CREE
            IF (ABS (X) .GE. EPS)
     :         CALL SPEAK (JA,JB,IA1,IB1,IB2,IA2,KK,X)
   40    CONTINUE
      ENDIF
*
   41 RETURN
*
  300 FORMAT('COR: Spectator quantum numbers not diagonal for',
     :       ' non-interacting shells')
  301 FORMAT('COR: Dimension error: reset PARAMETER IDIM to at least ',
     :       1I2,' and recompile.')
*
      END
