************************************************************************
*                                                                      *
      SUBROUTINE BREID (JA,JB,JA1,IPCA,JB1)
*                                                                      *
*   Computes closed shell contributions - aaaa and exchange only.      *
*                                                                      *
*   Call(s) to: [LIB92]: CLRX, CXK, ITRIG, TALK, SNRC.                 *
*                                                                      *
*                                           LAST UPDATE: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION CONE(7,20),JS(4),KAPS(4),KS(4),S(12)
*
      COMMON/BCORE/ICORE(NNNW)
     :      /CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /M1/NQ1(NNNW),NQ2(NNNW)
     :      /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB4/NP(NNNW),NAK(NNNW)
*
      PARAMETER (NUMAX = 20)
*
*   1.0  Initialization
*
      IF (IPCA .EQ. 2) THEN
         IA1 = KLIST(JA1)
      ELSE
         IA1 = JLIST(JA1)
      ENDIF
      IB1 = KLIST(JB1)
*
      ISG = 1
      IF (JA .EQ. JB) THEN
         IF ((ICORE(IA1) .NE. 0) .AND. (ICORE(IB1) .NE. 0)) THEN
            IF (JA .GT. 1) RETURN
            ISG = -1
         ENDIF
      ENDIF
*
      JS(1) = IA1
      JS(2) = IB1
      JS(3) = IA1
      JS(4) = IB1
      NQS1 = NQ1(IA1)
      NQS2 = NQ2(IB1)
      DO 1 I = 1,4
         KAPS(I) = 2*NAK(JS(I))
         KS(I) = IABS (KAPS(I))
    1 CONTINUE
      CONST = NQS1*NQS2
      IF (IBUG2 .NE. 0) WRITE (99,300) IA1,IB1
*
*   2.0  Set range of tensor indices
*
      CALL SNRC (JS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
      IF (IBUG2 .NE. 0) WRITE (99,301) ND1,ND2,NE1,NE2,IBRD,IBRE
      IF (IA1 .NE. IB1) GOTO 3
*
*   3.0 Calculate aaaa interaction
*
      DO 2 N = 1,ND2
         NU = ND1+2*(N-1)
         K = NU
         IF (MOD (K,2) .NE. 1) RETURN
         KAP1 = KAPS(1)/2
         GAM = CLRX (KAP1,NU,KAP1)
         DKSKS = KS(1)*KS(1)
         DNUNU1 = NU*(NU+1)
         COEF = CONST*TWO*DKSKS*GAM*GAM/DNUNU1
         IF (IBUG2 .NE. 0) WRITE (99,302) NU,GAM,COEF
         ITYPE = ISG*4
         CALL TALK (JA,JB,NU,IA1,IA1,IA1,IA1,ITYPE,COEF)
    2 CONTINUE
      RETURN
*
*   Calculate exchange interactions
*
    3 CONTINUE
      IF (IBRE .LT. 0) RETURN
      IF (NE2 .GT. NUMAX) THEN
         WRITE (*,304)
         STOP
      ENDIF
*
      DO 4 N = 1,NE2
         DO 4 MU = 1,7
            CONE(MU,N) = ZERO
    4 CONTINUE
*
      PROC = -CONST/DBLE (KS(1)*KS(2))
*
*   Negative sign arises from Pauli phase factor
*
      DO 10 N = 1,NE2
         NU = NE1+2*(N-1)
         K = NU
         IP = (KS(1)-KS(2))/2+K
         IPP = IP+1
         IF (NU .EQ. 0) GOTO 8
         KK = K+K+1
         IF (ITRIG(KS(1),KS(2),KK) .EQ. 0) GOTO 6
         PROD = PROC
         IF (MOD (IP,2) .NE. 0) PROD = -PROD
         CALL CXK (S,JS,KAPS,NU,K,IBRE,2)
         IF (IBUG2 .NE. 0) WRITE (99,303) PROD,(S(MU),MU = 1,3)
         DO 5 MU = 1,3
            CONE (MU,N) = CONE(MU,N)+PROD*S(MU)
    5    CONTINUE
*
    6    K = NU-1
         KK = K+K+1
         IF (ITRIG(KS(1),KS(2),KK) .EQ. 0) GOTO 8
         PROD = PROC
         IF (MOD (IPP,2) .NE. 0) PROD = -PROD
         CALL CXK (S,JS,KAPS,NU,K,IBRE,2)
         IF (IBUG2 .NE. 0) WRITE (99,303) PROD,(S(MU),MU = 1,3)
         DO 7 MU = 1,3
            CONE(MU,N) = CONE(MU,N)+PROD*S(MU)
    7    CONTINUE
*
    8    IF (N .EQ. NE2) GOTO 11
         K = NU+1
         KK = K+K+1
         PROD = PROC
         IF (MOD (IPP,2) .NE. 0) PROD = -PROD
         CALL CXK (S,JS,KAPS,NU,K,IBRE,2)
         IF (IBUG2 .NE. 0) WRITE (99,303) PROD,(S(MU),MU = 1,7)
         DO 9 MU = 1,7
            CONE (MU,N) = CONE(MU,N)+PROD*S(MU)
    9    CONTINUE
   10 CONTINUE
*
*   4.0  Output results
*
   11 CONTINUE
*
      DO 12 N = 1,NE2
         NU = NE1+2*(N-1)
         ITYPE = ISG*5
         CALL TALK (JA,JB,NU,IB1,IA1,IB1,IA1,ITYPE,CONE(1,N))
         CALL TALK (JA,JB,NU,IA1,IB1,IB1,IA1,ITYPE,CONE(2,N))
         CALL TALK (JA,JB,NU,IA1,IB1,IA1,IB1,ITYPE,CONE(3,N))
         IF (N .EQ. NE2) GOTO 12
         NUP1 = NU+1
         ITYPE = ISG*6
         CALL TALK (JA,JB,NUP1,IA1,IB1,IA1,IB1,ITYPE,CONE(4,N))
         CALL TALK (JA,JB,NUP1,IB1,IA1,IB1,IA1,ITYPE,CONE(5,N))
         CALL TALK (JA,JB,NUP1,IA1,IB1,IB1,IA1,ITYPE,CONE(6,N))
         CALL TALK (JA,JB,NUP1,IB1,IA1,IA1,IB1,ITYPE,CONE(7,N))
   12 CONTINUE
      RETURN
*
  300 FORMAT ('BREID: orbitals ',2I3)
  301 FORMAT (2X,'ND1 ND2 NE1 NE2 IBRD IBRE ',6I5)
  302 FORMAT (2X,'aaaa contribution: NU,GAM,COEF',I5,2(3X,1PD15.8))
  303 FORMAT (2X,'PROD = ',1PD15.8
     :       /' S',7D15.8)
  304 FORMAT ('BREID: Dimension error for NUMAX.')
*
      END
