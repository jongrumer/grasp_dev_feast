********************************************************************
*                                                                  *
      SUBROUTINE EL1(JJA,JJB,JA,JB,IIRE,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 03  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :           N'1 = N1        *
*                                                  N'2 = N2        *
*                                                                  *
*      SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,      *
*                         RECO,RECO2,SIXJ,SPEAK,WW1                *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
      COMMON /ORB4/NP(NNNW),NAK(NNNW)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      DIMENSION CONE(7,20),S(12),IS(4),KAPS(4),KS(4)
      DIMENSION PMGG(30),RAGG(30),J(2)
      IF(JA.NE.JB)GO TO 9
CGG
      IF(JJA.NE.JJB) RETURN
C
C     THE CASE 1111   + + - -
C
      IF(IIRE.EQ.0)GO TO 50
      CALL RECO(JA,JA,JA,JA,0,IAT)
      IF(IAT.EQ.0)RETURN
   50 CALL PERKO2(JA,JA,JA,JA,1)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      IA=JLIST(JA)
      J(1)=ID1(3)
      IP2=ITREXG(J(1),J(1),J(1),J(1),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      L1=(J(1)+1)/2
      IP1=IP2
      IG1=IG2
      IF (ICOLBREI .EQ. 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IA
        IS(4)=IA
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=ABS(KAPS(1))
        KS(2)=ABS(KAPS(2))
        KS(3)=ABS(KAPS(3))
        KS(4)=ABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD .LE. 0)RETURN 
      END IF
      DO 5 I2=IP1,IG1,2
        KRA=(I2-1)/2
        IF (ICOLBREI .EQ. 1) THEN
          CALL COULOM(L1,L1,L1,L1,ID1(5),ID1(5),ID1(5),ID1(5),KRA,A1)
          IF(ABS(A1).LT.EPS)GO TO 5
          A1=-A1*HALF
        END IF
      AB=ZERO
      DO 6 I3=IP2,IG2,2
        J12=(I3-1)/2
        IF(IXJTIK(J(1),J(1),KRA*2,J(1),J(1),J12*2).EQ.0)GO TO 6
        CALL WW1(IK1,BK1,ID1,BD1,J12,QM1,QM2,QM3,QM4,AA)
        IF(ABS(AA).LT.EPS)GO TO 6
        CALL SIXJ(J(1),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
        AA=AA*SI*SQRT(DBLE(I3))
        IFAZ=IK1(3)+J12+KRA
        IF((IFAZ/2)*2.NE.IFAZ)AA=-AA
        AB=AB+AA
    6 CONTINUE
C
C     RECOUPLING COEFFICIENTS
C
      IF (ICOLBREI .EQ. 1) THEN
        BB=AB*A1
        BB=BB/SQRT(DBLE(IK1(6)+1))
        IF(DABS(BB).GT.EPS)
     :           CALL SPEAK(JJA,JJB,IA,IA,IA,IA,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        N=(KRA-ND1)/2+1
        IF(((KRA-ND1)/2)*2 .EQ. (KRA-ND1)) THEN
          CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
          IF(DABS(S(1)).GT.EPS) THEN
            BB =-HALF*S(1)*AB/SQRT(DBLE(IK1(6)+1))
            IF(DABS(BB).GT.EPS) 
     :         CALL TALK(JJA,JJB,KRA,IA,IA,IA,IA,4,BB)
          END IF
        END IF
      END IF
    5 CONTINUE
      RETURN
C  ............................................................
    9 IF(NPEEL.LE.1)RETURN
      IF(IIRE.EQ.0)GO TO 51
      CALL RECO(JA,JB,JB,JB,1,IAT)
      IF(IAT.EQ.0)RETURN
   51 IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=-HALF
      QM3=HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      IP3=IP1
      IG3=IG1
      DO 11 I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1.GT.30)GO TO 10
        RAGG(KRA1)=ZERO
        PMGG(KRA1)=ZERO
        CALL RECO2(JA,JB,KRA*2,0,IAT,REC)
        IF(IAT.EQ.0)GO TO 11
        CALL GG1122(KRA,KRA,QM1,QM2,QM3,QM4,RAG)
        IF(ABS(RAG).LT.EPS) GO TO 11
        RAGG(KRA1)=RAG
        CALL RECO2(JA,JB,KRA*2,1,IAT,REC)
        PMGG(KRA1)=REC
   11 CONTINUE
C * * *                      * * *                      * * *
C     CASES 1212   + + - -        TRANSFORM TO  1122   + - + -
C           2121                                1122
C
      IF (ICOLBREI .EQ. 2) THEN
        IS(1)=IA
        IS(2)=IB
        IS(3)=IA
        IS(4)=IB
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=ABS(KAPS(1))
        KS(2)=ABS(KAPS(2))
        KS(3)=ABS(KAPS(3))
        KS(4)=ABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        DO 16 II=1,20
           CONE(1,II)=ZERO
           CONE(2,II)=ZERO
           CONE(3,II)=ZERO
           CONE(4,II)=ZERO
           CONE(5,II)=ZERO
           CONE(6,II)=ZERO
           CONE(7,II)=ZERO
   16   CONTINUE
        IF(IBRD .EQ. 0 .AND. IBRE .EQ.0)RETURN 
      END IF
      DO 1 I1=IP1,IG1,2
      KRA=(I1-1)/2
      KRA1=KRA+1
      IF(KRA1.GT.30)GO TO 10
      IF (ICOLBREI .EQ. 1) THEN
        CALL COULOM(L1,L2,L1,L2,ID1(5),ID2(5),ID1(5),ID2(5),KRA,AA)
        IF(ABS(AA).LT.EPS)GO TO 1
        AA=AA*PMGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 1
        AA=AA*RAGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 1
        AA=AA/SQRT(DBLE(I1))
        CALL SPEAK(JJA,JJB,IA,IB,IA,IB,KRA,AA)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        N=(KRA-ND1)/2+1
        IF(((KRA-ND1)/2)*2 .EQ. (KRA-ND1)) THEN
          CALL CXK(S,IS,KAPS,KRA,KRA,3,1)
          IF(DABS(S(1)).GT.EPS) THEN
            BB=S(1)*PMGG(KRA1)*RAGG(KRA1)/SQRT(DBLE(I1))
            IF(DABS(BB).GT.EPS)CALL TALK(JJA,JJB,KRA,IA,IA,IB,IB,4,BB)
          END IF
        END IF
      END IF
    1 CONTINUE
C * * *                      * * *                      * * *
C     CASES 1221   + + - -        TRANSFORM TO  1122   + - + -
C           2112                                1122
C
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      DO 2 I2=IP2,IG2,2
      KRA=(I2-1)/2
      IF(KRA.GT.30)GO TO 10
      IF (ICOLBREI .EQ. 1) THEN
        CALL COULOM(L1,L2,L2,L1,ID1(5),ID2(5),ID2(5),ID1(5),KRA,A1)
        IF(ABS(A1).LT.EPS)GO TO 2
      END IF
      AB=ZERO
      DO 3 I3=IP3,IG3,2
        J12=(I3-1)/2
        KRA1=J12+1
        IF(KRA1.GT.30)GO TO 10
        AA=PMGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 3
        AA=AA*RAGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 3
        IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2).EQ.0)GO TO 3
        CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
        AA=AA*SI*SQRT(DBLE(I3))
        AB=AB+AA
    3 CONTINUE
      IF (ICOLBREI .EQ. 1) THEN
        BB=A1*AB
        IF(ABS(BB).GT.EPS)CALL SPEAK(JJA,JJB,IA,IB,IB,IA,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        NU=KRA 
        IF(((NU-NE1)/2)*2 .EQ. (NU-NE1)) THEN
          IF((ITRIG(KS(1),KS(4),NU+NU+1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(3),NU+NU+1).NE.0)) THEN
            IF(NU .GT. 0) THEN
              N=(NU-NE1)/2+1
              CALL CXK(S,IS,KAPS,NU,KRA,4,2)
              DO 21 MU = 1,3
                CONE(MU,N)=CONE(MU,N)+AB*S(MU)
   21         CONTINUE
            END IF
          END IF
        END IF
        NU=KRA+1
        IF(((NU-NE1)/2)*2 .EQ. (NU-NE1)) THEN
          IF((ITRIG(KS(1),KS(4),NU+NU-1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(3),NU+NU-1).NE.0)) THEN
            IF(NU .GE. 0) THEN
              N=(NU-NE1)/2+1
              IF(N .LE. NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                DO 22 MU = 1,3
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
   22           CONTINUE
              END IF
            END IF
          END IF
        END IF
        NU=KRA-1
        IF(((NU-NE1)/2)*2 .EQ. (NU-NE1)) THEN
          IF((ITRIG(KS(1),KS(4),NU+NU+3).NE.0) .AND.
     :       (ITRIG(KS(2),KS(3),NU+NU+3).NE.0)) THEN
            IF(NU .GE. 0) THEN
              N=(NU-NE1)/2+1
              IF(N .LT. NE2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,4,2)
                DO 23 MU = 1,7
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
   23           CONTINUE
              END IF
            END IF
          END IF
        END IF
      END IF
    2 CONTINUE
      IF (ICOLBREI .EQ. 2) THEN
        DO 36 N = 1,NE2
           NU=NE1+2*(N-1)
           CALL TALK(JJA,JJB,NU,IB,IA,IB,IA,5,CONE(1,N))
           CALL TALK(JJA,JJB,NU,IA,IB,IB,IA,5,CONE(2,N))
           CALL TALK(JJA,JJB,NU,IA,IB,IA,IB,5,CONE(3,N))
           IF(N.EQ.NE2) GO TO 36
           NUP1=NU+1
           CALL TALK(JJA,JJB,NUP1,IA,IB,IA,IB,6,CONE(4,N))
           CALL TALK(JJA,JJB,NUP1,IB,IA,IB,IA,6,CONE(5,N))
           CALL TALK(JJA,JJB,NUP1,IA,IB,IB,IA,6,CONE(6,N))
           CALL TALK(JJA,JJB,NUP1,IB,IA,IA,IB,6,CONE(7,N))
   36   CONTINUE
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL1  PMGG RAGG')
      STOP
      END
