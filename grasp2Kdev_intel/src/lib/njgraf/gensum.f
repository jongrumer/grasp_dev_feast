************************************************************************
*                                                                      *
      SUBROUTINE GENSUM (J6C,J7C,J8C,J9C,JWC,J6,J7,J8,J9,JW,JDEL,
     :                   LDEL,SUMVAR,MP,J6P,J7P,J8P,J9P,JWORD,NLSUM,
     :                   NBJ,NB6J,K6CP,K7CP,K8CP,K9CP,JSUM4,JSUM5,
     :                   JSUM6,INV6J,RECUP)
*                                                                      *
*   Carries  out the summation over coefficients defined by the arr-   *
*   ays J6, J7, J8, LDEL and  JW  to give RECUP. The entry is either   *
*   made from NJGRAF or directly  assuming that the arrays J6,...,JW   *
*   have already been determined  by  a previous entry to NJGRAf and   *
*   that the summation is required for another set of j values defi-   *
*   ned by the array J1. RECUP is the recoupling coefficient.          *
*                                                                      *
*   Call(s) to: [NJGRAF]: DRACAH, RDIAG.                               *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      LOGICAL LDIAG,NOEL,FREE,SUMVAR
*
      PARAMETER (
     :   MANGM = 60,  M3MNGM = 3*MANGM, MANGMP = 2*(MANGM/3),
     :   MTRIAD = 12, M2TRD = 2*MTRIAD, M4TRD = 4*MTRIAD,
     :   M6J = 20,MSUM = 10,MFACT = 500)
*
      PARAMETER (EPSIL = 1.0D-10)
c      INTEGER*8 IJ6CP
*
      DIMENSION MAT(MTRIAD,MTRIAD),NOEL(MTRIAD),
     :   MAXLP(MTRIAD),JSUM2(MTRIAD),JSUM3(MTRIAD),
     :   JSUM(2,M6J),JWTEST(M6J),WSTOR(M6J),IPAIR(2,2),LDIAG(MTRIAD)
      DIMENSION XJ1(MANGM),IST(6)
      DIMENSION J12(4,MTRIAD,MTRIAD)
      DIMENSION J6P(MANGMP),J7P(MANGMP),J8P(MANGMP),J9P(MANGMP),
     :   JWORD(6,M6J),NBJ(MSUM),NB6J(MSUM),K6CP(MSUM),
     :   K7CP(MSUM),K8CP(MSUM),K9CP(MSUM),JSUM6(MTRIAD),
     :   JSUM4(MTRIAD,M6J),JSUM5(MTRIAD,M6J),INV6J(M6J)
      DIMENSION J6(M3MNGM),J7(M3MNGM),J8(M3MNGM),
     :   J9(MANGMP),JW(6,M6J),LDEL(M6J,2),SUMVAR(MANGM)
*
      COMMON/DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /FACTS/GAM(MFACT)
     :      /COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),FREE(MANGM)
*
      DATA MXCSVR/4/
*
*   evaluates all terms in J6, J7, J8, J9, LDEL, JW which do not
*   involve a summation. The result is stored in RECUP and IASTOR
*
      IF (IBUG3 .EQ. 1) THEN
*
         DO 139 I = 1,M
            XJ1(I) = 0.5D 00*DBLE (J1(I)-1)
  139    CONTINUE
*
         WRITE (99,400) (XJ1(I),I = 1,M)
         WRITE (99,306) NLSUM
         WRITE (99,401)
      ENDIF
*
      MM = M+1
      J1(MM) = 1
*
*   Test delta functions
*
      J1(MM) = 1
      IF (JDEL .LE. 0) GOTO 180
*
      DO 141 I = 1,JDEL
         I1 = LDEL(I,1)
         I2 = LDEL(I,2)
         IF ((I1 .GT. MM) .OR. (I2 .GT. MM))THEN
            IF (I1.GT.MM) J1(I1) = J1(I2)
            IF (I2.GT.MM) J1(I2) = J1(I1)
         ELSE
            IF (J1(I1) .NE. J1(I2)) THEN
               RECUP = 0.0D 00
               RETURN
            ENDIF
         ENDIF
  141 CONTINUE
*
  180 RECUP = 1.0D 00
      IF (JWC .NE. 0) THEN
*
*   Multiply RECUP by all Racah coefficients which do not involve a
*   summation
*
         IF (IBUG3 .EQ. 1) WRITE (99,309)
*
         DO 7 I = 1,JWC
            IF (INV6J(I) .GT. 0) GOTO 7
            DO 3 J = 1,6
               I1 = JW(J,I)
               IST(J) = J1(I1) - 1
    3       CONTINUE
*
            CALL DRACAH (IST(1),IST(2),IST(3),IST(4),IST(5),IST(6),X1)
            IF (IBUG3 .EQ. 1) WRITE (99,305) (XJ1(JW(K,I)),K = 1,6),X1
            RECUP = RECUP*X1
*
    7    CONTINUE
*
      ENDIF
*
      SQR = 1.0D 00
*
      IF (J6C .NE. 0) THEN
         DO 12 I = 1,J6C
            I1 = J6(I)
            SQR = SQR*J1(I1)
   12     CONTINUE
      ENDIF
*
      SPR = 1.0D 00
*
      IF (J9C .NE. 0) THEN
         DO 144 I = 1,J9C
            I1 = J9(I)
            SPR = SPR*J1(I1)
  144    CONTINUE
      ENDIF
*
      RECUP = RECUP*SQRT (SQR/SPR)
      IF (ABS(RECUP) .LT. EPSIL) GOTO 145
      IASTOR = 0
*
      IF (J7C .NE. 0) THEN
         DO 17 I = 1,J7C
            I1 = J7(I)
            IASTOR = IASTOR + J1(I1) -1
   17    CONTINUE
      ENDIF
*
      IF (J8C .NE. 0) THEN
         DO 22 I = 1,J8C
            I1 = J8(I)
            IASTOR = IASTOR +2*(J1(I1)-1)
   22    CONTINUE
      ENDIF
*
      IF (NLSUM .LE. 0) THEN
         IASTOR = IASTOR/2
*
*   No summation involved. End of computation
*
         STOR1 = 1.0D 00
         STOR = 1.0D 00
         IF (MOD (IASTOR,2) .EQ. 1) RECUP = -RECUP
         IF (IBUG3 .EQ. 1) WRITE (99,303) RECUP
         RETURN
*
      ENDIF
*
*   Evaluation of the part involving summations.
*
      NFS = 0
      JWR = 0
      J6F = 0
      J7F = 0
      J8F = 0
      J9F = 0
      NPS = 0
   25 NPS = NPS+1
      IF (IBUG3 .EQ. 1) WRITE (99,302) NPS
*
*   Loop on the disconnected summations
*
      IAS = 0
      NSUM = NBJ(NPS)-NFS
      JWRD = NB6J(NPS)-JWR
      J6CP = K6CP(NPS)
      J7CP = K7CP(NPS)
      J8CP = K8CP(NPS)
      J9CP = K9CP(NPS)
*
*   The range of values of each summation variable is defined by
*   establishing a matrix of the links between variables.
*   MAT(I,J) contains:
*       I = J  Number of possible values of I due to triangular
*              relations with non-variables, i.e. constants.
*       I > J  Number of links between I and J through constants
*       I < J  Value of the constant, if the above is 1. If not,
*              these values are srored in J12(L,I,J) where there
*              is room for MXCSVR such values (L .LE. 4)
*
      DO 52 I = 1,NSUM
         DO 152 J = 1,NSUM
            MAT(I,J) = 0
  152    CONTINUE
   52 CONTINUE
*
      DO 66 I1 = 1,NSUM
         I1T = I1+NFS
         I2 = JSUM6(I1T)
         DO 65 I3 = 1,I2
            I = JSUM5(I1T,I3)
            J = JSUM4(I1T,I3)
            GOTO (54,55,56,57,58,59),J
*
*   The rows of the IPAIR arrays give limits of summation imposed
*
   54       IPAIR(1,1) = JWORD(2,I)
            IPAIR(1,2) = JWORD(5,I)
            IPAIR(2,1) = JWORD(3,I)
            IPAIR(2,2) = JWORD(6,I)
            GOTO 60
*
   55       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(5,I)
            IPAIR(2,1) = JWORD(4,I)
            IPAIR(2,2) = JWORD(6,I)
            GOTO 60
*
   56       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(6,I)
            IPAIR(2,1) = JWORD(4,I)
            IPAIR(2,2) = JWORD(5,I)
            GOTO 60
*
   57       IPAIR(1,1) = JWORD(2,I)
            IPAIR(1,2) = JWORD(6,I)
            IPAIR(2,1) = JWORD(3,I)
            IPAIR(2,2) = JWORD(5,I)
            GOTO 60
*
   58       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(2,I)
            IPAIR(2,1) = JWORD(3,I)
            IPAIR(2,2) = JWORD(4,I)
            GOTO 60
*
   59       IPAIR(1,1) = JWORD(1,I)
            IPAIR(1,2) = JWORD(3,I)
            IPAIR(2,1) = JWORD(2,I)
            IPAIR(2,2) = JWORD(4,I)
*
   60       DO 63 I4 = 1,2
               KM = 0
               DO 62 I5 = 1,2
                  IF (IPAIR(I4,I5) .GT. MP) KM = KM+1
   62          CONTINUE
*
               JJ1 = IPAIR(I4,1)
               JJ2 = IPAIR(I4,2)
               IF (KM .EQ. 1) GOTO 67
               IF (KM .GT. 1) GOTO 63
*
*   One variable linked to two constants. Fix the diagonal MAT(I,I)
*
              JT1 = J1(JJ1)-1
              JT2 = J1(JJ2)-1
              JMIN = ABS (JT1-JT2)
              JMAX = JT1+JT2
*
              IF (MAT(I1,I1) .GT. 1) THEN
*
*   If there are several couples of constants, take the more
*   stringent combination
*
                 JMIN = MAX (JMIN,JSUM(1,I1))
                 JMAX = MIN (JMAX,JSUM(2,I1))
                 IF (JMAX .GE. JMIN) THEN
                    JSUM(1,I1) = JMIN
                    JSUM(2,I1) = JMAX
                    MAT(I1,I1) = (JMAX-JMIN)/2+1
                    GOTO 63
                 ELSE
                    RECUP = 0.0D 00
                    GOTO 110
                 ENDIF
              ELSEIF (MAT(I1,I1) .LT. 1) THEN
*
*   First time
*
                  MAT(I1,I1) = (JMAX-JMIN)/2+1
                  JSUM(1,I1) = JMIN
                  JSUM(2,I1) = JMAX
               ENDIF
*
               GOTO 63
*
*   One variable linked to one constant and one variable  non diagonal
*   element
*
   67          JT1 = MIN (JJ1,JJ2)
               JT2 = MAX (JJ1,JJ2)-MP
               IF (JT2 .GT. I1) GOTO 63
               JT4 = J1(JT1)-1
               K = MAT(I1,JT2)
               IF (K .EQ. 0) GOTO 107
*
               DO 71 LL = 1,K
                  IF (JT4 .EQ. J12(LL,JT2,I1)) GOTO 63
   71          CONTINUE
*
  107          K = K+1
               IF (K .GT. MXCSVR) GOTO 63
               MAT(I1,JT2) = K
               J12(K,JT2,I1) = JT4
*
   63       CONTINUE
   65    CONTINUE
   66 CONTINUE
*
*   Reduce the diagonal elements by taking into account the non
*   diagonal elements, and keep the latter only if needed
*
  150 ICHAN = 0
*
      DO 74 I = 1,NSUM
         NOEL(I) = .TRUE.
         I1 = I-1
         IF (I1 .EQ. 0) GOTO 170
         DO 72  J = 1,I1
            IF ((MAT(I,J) .EQ. 0) .OR. (MAT(J,J) .EQ. 0)) GOTO 72
            IK1 = I
            IK2 = J
            CALL RDIAG (I,J,IK1,IK2,ICHAN,MAT,JSUM,J12)
            NOEL(I) = .FALSE.
   72    CONTINUE
  170    IF (I .EQ. NSUM) GOTO 74
         I2 = I+1
*
         DO 73 J = I2,NSUM
            IF ((MAT(J,I) .EQ. 0) .OR. (MAT(J,J) .EQ. 0)) GOTO 73
            IK1 = J
            IK2 = I
            CALL RDIAG (I,J,IK1,IK2,ICHAN,MAT,JSUM,J12)
   73    CONTINUE
   74 CONTINUE
*
      IF (ICHAN .NE. 0) GOTO 150
      GOTO 220
*
*   Carry out the summations.
*
  220 DO 230 I = 1,NSUM
         JSUM3(I) = 1
         LDIAG(I) = .FALSE.
         IF (MAT(I,I) .EQ. 1) LDIAG(I) = .TRUE.
  230 CONTINUE
*
      DO 231 I = 1,JWRD
         JWTEST(I) = 1
  231 CONTINUE
*
      STOR = 0.0D 00
      STOR1 = 1.0D 00
      NOLP = 0
      IP = 1
  240 NOLP = NOLP+1
*
*   Find the range of JSUM2(NOLP)
*   NOLP is the index  of the summation variable
*
      JMIN = JSUM(1,NOLP)
      JMAX = JSUM(2,NOLP)
      IF (NOEL(NOLP)) GOTO 241
      NO1 = NOLP-1
*
      DO 242 NJ = 1,NO1
         IF (MAT(NOLP,NJ) .EQ. 1) THEN
            JJ1 = MAT(NJ,NOLP)
            JJ2 = JSUM2(NJ)
            JMIN = MAX (JMIN,IABS(JJ2-JJ1))
            JMAX = MIN (JMAX,JJ1+JJ2)
         ELSEIF (MAT(NOLP,NJ) .GT. 1) THEN
            K = MAT(NOLP,NJ)
            JJ2 = JSUM2(NJ)
*
            DO 245 I = 1,K
            JJ1 = J12(I,NJ,NOLP)
            JMIN = MAX (JMIN,IABS(JJ2-JJ1))
            JMAX = MIN (JMAX,JJ1+JJ2)
  245       CONTINUE
*
         ENDIF
*
  242 CONTINUE
*
  241 JSUM2(NOLP) = JMIN
      MAXLP(NOLP) = JMAX
      IF (LDIAG(NOLP)) JSUM3(NOLP) = 0
      IF (NOLP .LT. NSUM) GOTO 240
*
      DO 260 JJ = JMIN,JMAX,2
         JSUM2(NSUM) = JJ
*
*   Determine which RACAH coefficients need re-evaluating and
*   set JWTEST appropriately
*
      DO 114 J = IP,NSUM
         IF (JSUM3(J) .LE. 0) GOTO 114
         I2 = JSUM6(J)
*
         DO 113 I1 = 1,I2
            I3 = JSUM5(J,I1)
            JWTEST(I3) = 1
  113    CONTINUE
  114 CONTINUE
*
      DO 98 J = 1,JWRD
         IF (JWTEST(J) .EQ. 0) GOTO 98
         JWJ = J+JWR
*
         DO 90 I = 1,6
            IF (JWORD(I,JWJ) .LE. MP) THEN
               I1 = JWORD(I,JWJ)
               IST(I) = J1(I1) - 1
            ELSE
               I1 = JWORD(I,JWJ)-MP-NFS
               IST(I) = JSUM2(I1)
            ENDIF
   90    CONTINUE
*
         CALL DRACAH (IST(1),IST(2),IST(3),IST(4),IST(5),IST(6),X1)
         WSTOR(J) = X1
         IF (IBUG3 .EQ. 1) THEN
            DO 99 I = 1,6
               XJ1(I) = 0.5D 00*DBLE (IST(I))
   99       CONTINUE
*
            WRITE (99,305) (XJ1(I), I = 1,6),X1
         ENDIF
   98 CONTINUE
*
*   Form product of Racah coefficients, (2J+1) factors and (-1)
*   factors in STOR1
*
      DO 126 I = 1,JWRD
         STOR1 = STOR1*WSTOR(I)
  126 CONTINUE
*
*   IASTOR contains the power of (-1) which is common to all terms
*
      IX2 = 0
c     IJ6CP = 1
      DIJ6CP = 1.0
      IF (J6CP .NE. J6F) THEN
         JB = J6F+1
*
         DO 128 I = JB,J6CP
            I1 = J6P(I)-NFS
C      print *,IJ6CP,JSUM2(I1),i1,'zou,gensum,1'
c           IJ6CP = IJ6CP*(JSUM2(I1)+1)
            DIJ6CP = DIJ6CP*DBLE(JSUM2(I1)+1)
  128    CONTINUE
      ENDIF
*
      IF (J9CP .NE. J9F) THEN
         JB = J9F+1
*
         DO 147 I = JB,J9CP
            I1 = J9P(I)-NFS
C      print *,IJ6CP,JSUM2(I1),i1,'zou,gensum,2'
c           IJ6CP = IJ6CP/(JSUM2(I1)+1)
            DIJ6CP = DIJ6CP/DBLE(JSUM2(I1)+1)
  147    CONTINUE
      ENDIF
*
C      print *,STOR1, IJ6CP,'zou,gensum'
c     STOR1 = STOR1*SQRT (DBLE (IJ6CP))
      STOR1 = STOR1*SQRT (DIJ6CP)
*
      IF (J7CP .NE. J7F) THEN
         JB = J7F+1
*
         DO 131 I = JB,J7CP
            I1 = J7P(I)-NFS
            IX2 = IX2 + JSUM2(I1)
  131    CONTINUE
      ENDIF
*
      IF (J8CP .NE. J8F) THEN
         JB = J8F+1
*
         DO 134 I = JB,J8CP
            I1 = J8P(I)-NFS
            IX2 = IX2 + 2*(JSUM2(I1))
  134    CONTINUE
      ENDIF
*
      IF (MOD(IX2,2) .EQ. 1) THEN
         IAS = -1
         IX2 = IX2+1
      ENDIF
*
      IX2 = IX2/2
*
*   Add term into STOR and reset STOR1 to 1 ready for next term
*
      IF (MOD(IX2,2) .EQ. 1) STOR1 = -STOR1
      STOR = STOR + STOR1
      STOR1 = 1.0D 00
      NSUM1 = NSUM-1
      IF (NSUM1 .EQ. 0) GOTO 260
*
      DO 261 IK = 1,NSUM1
         JSUM3(IK) = 0
  261 CONTINUE
*
      DO 262 IK = 1,JWRD
         JWTEST(IK) = 0
  262 CONTINUE
*
  260 CONTINUE
*
  250 NOLP = NOLP-1
*
      IF (NOLP .NE. 0) THEN
         IF (LDIAG(NOLP)) GOTO 250
         JSUM3(NOLP) = 1
         JSUM2(NOLP) = JSUM2(NOLP)+2
         IF (JSUM2(NOLP) .GT. MAXLP(NOLP)) GOTO 250
         IP = NOLP
*
*   Proceed to next variable
*
         GOTO 240
*
      ENDIF
*
      RECUP = RECUP*STOR
      IF (IBUG3 .EQ. 1) WRITE (99,307) NPS,STOR,RECUP
      IF (ABS(RECUP) .LT. EPSIL) GOTO 145
      JWR = JWRD+JWR
      NFS = NSUM+NFS
      J6F = J6CP
      J7F = J7CP
      J8F = J8CP
      J9F = J9CP
      IASTOR = IASTOR+IAS
*
*   Proceed to next sum
*
      IF (NPS .LT. NLSUM) GOTO 25
      IASTOR = IASTOR/2
      IF (MOD (IASTOR,2) .NE. 0) RECUP = -RECUP
      IF (IBUG3 .EQ. 1) WRITE (99,304) RECUP
  110 RETURN
*
*   No summations. Check that there are no inconsistencies. Then
*   multiply by (-1) factor and exit
*
  145 RECUP = 0.0D 00
      RETURN
*
  302 FORMAT (' Sum Nr.',I3)
  303 FORMAT (' No summation. Recoupling coefficient = ',G15.8)
  304 FORMAT (' Recoupling coefficient = ',G15.8)
  305 FORMAT (6F5.1,10X,G15.8)
  306 FORMAT (' Number of independent sums:',I3)
  307 FORMAT (' Sum Nr.',I2,' Sum value = ',G15.8,' RECUP = ',G15.8)
  309 FORMAT (' Not involving summation variable')
  400 FORMAT (//' Printout from SUBROUTINE GENSUM'
     :       //' Values of angular momenta in *REAL* FORMAT'
     :        /(14F5.1))
  401 FORMAT (/' Racah W functions (6J)'
     :       /' Arguments in *REAL* FORMAT',18X,'value')
*
      END
