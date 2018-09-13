********************************************************************
*                                                                  *
      SUBROUTINE EL53(JJA,JJB,JA,JB,JC,JD,IREZ,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 14  -------------  *
*                                                                  *
*     THIS PACKAGE EVALUATED THE CASES - 1423, 4132, 1432, 4123    *
*                                                   ( IREZ = 1),   *
*                                                   ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 - 1   *
*                                                   N'2 = N2 + 1   *
*                                                   N'3 = N3 + 1   *
*                                                   N'4 = N4 - 1   *
*     AND    2314, 3241, 2341, 3214                 ( IREZ = 2),   *
*                                                   ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 + 1   *
*                                                   N'2 = N2 - 1   *
*                                                   N'3 = N3 - 1   *
*                                                   N'4 = N4 + 1   *
*                                                                  *
*     SUBROUTINE CALLED: COULOM,GG1234,ITREXG,IXJTIK,PERKO2,       *
*                      RECO,REC4,SIXJ,SPEAK                        *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/M1/NQ1(NNNW),NQ2(NNNW)
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
      COMMON /ORB4/NP(NNNW),NAK(NNNW)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      DIMENSION PMGG(30),J(4)
      DIMENSION COND(12,20),CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
      IF(NPEEL.LE.3)RETURN
      CALL RECO(JA,JD,JC,JB,3,IAT)
      IF(IAT.EQ.0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      ID=JLIST(JD)
      IF(IREZ.EQ.2)GO TO 20
      QM1=HALF
      QM2=-HALF
      QM3=-HALF
      QM4=HALF
      GO TO 21
   20 QM1=-HALF
      QM2=HALF
      QM3=HALF
      QM4=-HALF
   21 CALL PERKO2(JA,JB,JC,JD,4)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      J(4)=ID4(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      L4=(J(4)+1)/2
      IF (ICOLBREI .EQ. 2) THEN
        IF(IREZ .EQ. 1) THEN
          IS(1)=IA
          IS(2)=ID
          IS(3)=IB
          IS(4)=IC
        ELSE IF(IREZ .EQ. 2) THEN
          IS(1)=IB
          IS(2)=IC
          IS(3)=IA
          IS(4)=ID
        END IF
        KAPS(1)=2*NAK(IS(1))
        KAPS(2)=2*NAK(IS(2))
        KAPS(3)=2*NAK(IS(3))
        KAPS(4)=2*NAK(IS(4))
        KS(1)=ABS(KAPS(1))
        KS(2)=ABS(KAPS(2))
        KS(3)=ABS(KAPS(3))
        KS(4)=ABS(KAPS(4))
        CALL SNRC(IS,KAPS,KS,ND1,ND2,NE1,NE2,IBRD,IBRE)
        IF(IBRD .LE. 0 .AND. IBRE .LE. 0)RETURN
        DO 8 II=1,20
          COND(1,II) =ZERO
          COND(2,II) =ZERO
          COND(3,II) =ZERO
          COND(4,II) =ZERO
          COND(5,II) =ZERO
          COND(6,II) =ZERO
          COND(7,II) =ZERO
          COND(8,II) =ZERO
          COND(9,II) =ZERO
          COND(10,II)=ZERO
          COND(11,II)=ZERO
          COND(12,II)=ZERO
          CONE(1,II) =ZERO
          CONE(2,II) =ZERO
          CONE(3,II) =ZERO
          CONE(4,II) =ZERO
          CONE(5,II) =ZERO
          CONE(6,II) =ZERO
          CONE(7,II) =ZERO
          CONE(8,II) =ZERO
          CONE(9,II) =ZERO
          CONE(10,II)=ZERO
          CONE(11,II)=ZERO
          CONE(12,II)=ZERO
    8   CONTINUE
      END IF
      CALL GG1234(IK1,IK2,IK3,IK4,BK1,BK2,BK3,BK4,ID1,ID2,
     *ID3,ID4,BD1,BD2,BD3,BD4,QM1,QM2,QM3,QM4,RAG)
      IF(ABS(RAG).LT.EPS) RETURN
      IP1=ITREXG(J(1),J(2),J(3),J(4),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      DO 4 I4=IP1,IG1,2
      KRA=(I4-1)/2
      KRA1=KRA+1
      IF(KRA1.GT.30)GO TO 10
      PMGG(KRA1)=ZERO
      CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,0,
     *IAT,REC)
      IF(IAT.EQ.0)GO TO 4
      CALL RECO4(JA,JB,JC,JD,J(1),J(2),J(3),J(4),KRA*2,1,
     *IAT,REC)
      PMGG(KRA1)=REC
    4 CONTINUE
C
C     TRANSFORM FANO & RACAH PHASE CONVENTION
C     TO CONDON & SHORTLEY PHASE CONVENTION
C
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
     *+IK3(5)*IK3(4)-ID3(5)*ID3(4)+IK4(5)*IK4(4)-ID4(5)*ID4(4)
      IF((IFAZ/4)*4.NE.IFAZ)IFAZFRCS=-IFAZFRCS
C
      NN=0
      JB1=JB-1
      IFAZP=1
      DO 16 II=JA,JB1
      IN=JLIST(II)
      NN=NQ1(IN)+NN
   16 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
      NN=0
      JD1=JD-1
      DO 17 II=JC,JD1
      IN=JLIST(II)
      NN=NQ1(IN)+NN
   17 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
C * * *                      * * *                      * * *
C     CASES 1423   + + - -        TRANSFORM TO  1234   + - - +
C           4132                                1234
C                                                    (IREZ = 1)
C     OR
C     CASES 2314   + + - -        TRANSFORM TO  1234   - + + -
C           3241                                1234
C                                                    (IREZ = 2)
      DO 2 I3=IP1,IG1,2
      KRA=(I3-1)/2
      IF (ICOLBREI .EQ. 1) THEN
        IF(IREZ.EQ.2)GO TO 22
        CALL COULOM(L1,L4,L2,L3,ID1(5),ID4(5),ID2(5),ID3(5),KRA,A1)
        GO TO 23
   22   CALL COULOM(L2,L3,L1,L4,ID2(5),ID3(5),ID1(5),ID4(5),KRA,A1)
   23 IF(ABS(A1).LT.EPS)GO TO 2
      END IF
      KRA1=KRA+1
      IF(KRA1.GT.30)GO TO 10
      AA=PMGG(KRA1)
      IF(ABS(AA).LT.EPS) GO TO 2
      AA=AA*RAG
      IF(ABS(AA).LT.EPS) GO TO 2
      AA=AA/SQRT(DBLE(I3))
      IFAZ=J(4)+J(3)-2*KRA+2
      IF(IREZ.EQ.2)IFAZ=J(1)+J(2)-2*KRA+2
      IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
      AB=AA*DBLE(IFAZP)
      IF (ICOLBREI .EQ. 1) THEN
        BB=A1*AB*DBLE(IFAZFRCS)
        IF(IREZ.EQ.1)CALL SPEAK(JJA,JJB,IA,ID,IB,IC,KRA,BB)
        IF(IREZ.EQ.2)CALL SPEAK(JJA,JJB,IB,IC,IA,ID,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        NU=KRA
        IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
          IF((ITRIG(KS(1),KS(3),NU+NU+1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(4),NU+NU+1).NE.0)) THEN
            N=(NU-ND1)/2+1
            IF(NU .GT. 0) THEN
              CALL CXK(S,IS,KAPS,NU,KRA,1,1)
              DO 31 MU = 1,4
                COND(MU,N)=COND(MU,N)+AB*S(MU)
   31         CONTINUE
            END IF
          END IF
        END IF
        NU=KRA+1
        IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
          IF((ITRIG(KS(1),KS(3),NU+NU-1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(4),NU+NU-1).NE.0)) THEN
            N=(NU-ND1)/2+1
            IF(N .LE. ND2) THEN
              CALL CXK(S,IS,KAPS,NU,KRA,1,1)
              DO 32 MU = 1,4
                COND(MU,N)=COND(MU,N)+AB*S(MU)
   32         CONTINUE
            END IF
          END IF
        END IF
        NU=KRA-1
        IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
          IF((ITRIG(KS(1),KS(3),NU+NU+3).NE.0) .AND.
     :       (ITRIG(KS(2),KS(4),NU+NU+3).NE.0)) THEN
            IF(NU .GE. 0) THEN
              N=(NU-ND1)/2+1
              IF(N .LT. ND2) THEN
                CALL CXK(S,IS,KAPS,NU,KRA,1,1)
                DO 33 MU = 1,12
                  COND(MU,N)=COND(MU,N)+AB*S(MU)
   33           CONTINUE
              END IF
            END IF
          END IF
        END IF
      END IF
    2 CONTINUE
      IF (ICOLBREI .EQ. 2) THEN
        DO 36 N = 1,ND2
          NU=ND1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJA,JJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N.EQ.ND2) GO TO 36
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
   36   CONTINUE
      END IF
C * * *                      * * *                      * * *
C     CASES 1432   + + - -        TRANSFORM TO  1234   + - - +
C           4132                                1234
C                                                    (IREZ = 1)
C     OR
C     CASES 2341   + + - -        TRANSFORM TO  1234   - + + -
C           3214                                1234
C                                                    (IREZ = 2)
      IP2=ITREXG(J(1),J(3),J(4),J(2),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      DO 12 I2=IP2,IG2,2
      KRA=(I2-1)/2
      IF (ICOLBREI .EQ. 1) THEN
        IF(IREZ.EQ.2)GO TO 25
        CALL COULOM(L1,L4,L3,L2,ID1(5),ID4(5),ID3(5),ID2(5),KRA,A1)
        GO TO 24
   25   CALL COULOM(L2,L3,L4,L1,ID2(5),ID3(5),ID4(5),ID1(5),KRA,A1)
   24   IF(ABS(A1).LT.EPS)GO TO 12
      END IF
      AB=ZERO
      DO 13 I3=IP1,IG1,2
      J12=(I3-1)/2
      IFAZ=J(1)+J(4)+2*J12
      IF(IREZ.EQ.2)IFAZ=J(2)+J(3)+2*J12
      IF((IFAZ/2)*2.NE.IFAZ)GO TO 13
      KRA1=J12+1
      IF(KRA1.GT.30)GO TO 10
      AA=PMGG(KRA1)
      IF(ABS(AA).LT.EPS) GO TO 13
      AA=AA*RAG
      IF(ABS(AA).LT.EPS) GO TO 13
      IF(IXJTIK(J(1),J(3),KRA*2,J(4),J(2),J12*2).EQ.0)GO TO 13
      CALL SIXJ(J(1),J(3),KRA*2,J(4),J(2),J12*2,0,SI)
      AA=AA*SI*SQRT(DBLE(I3))
      IFAZ=J(3)-J(4)-4*KRA+2*J12
      IF(IREZ.EQ.2)IFAZ=J(1)-J(2)-4*KRA-2*J12
      IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
      AB=AB+AA
   13 CONTINUE
      IF(ABS(AB).LT.EPS)GO TO 12
      AB=AB*DBLE(IFAZP)
      IF (ICOLBREI .EQ. 1) THEN
        BB=A1*AB*DBLE(IFAZFRCS)
        IF(IREZ.EQ.1)CALL SPEAK(JJA,JJB,IA,ID,IC,IB,KRA,BB)
        IF(IREZ.EQ.2)CALL SPEAK(JJA,JJB,IB,IC,ID,IA,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        NU=KRA
        IF(((NU-NE1)/2)*2 .EQ. (NU-NE1)) THEN
          IF((ITRIG(KS(1),KS(4),NU+NU+1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(3),NU+NU+1).NE.0)) THEN
            N=(NU-NE1)/2+1
            IF(NU .GT. 0) THEN
              CALL CXK(S,IS,KAPS,NU,KRA,1,2)
              DO 41 MU = 1,4
                CONE(MU,N)=CONE(MU,N)+AB*S(MU)
   41         CONTINUE
            END IF
          END IF
        END IF
        NU=KRA+1
        IF(((NU-NE1)/2)*2 .EQ. (NU-NE1)) THEN
          IF((ITRIG(KS(1),KS(4),NU+NU-1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(3),NU+NU-1).NE.0)) THEN
            N=(NU-NE1)/2+1
            IF(N .LE. NE2) THEN
              CALL CXK(S,IS,KAPS,NU,KRA,1,2)
              DO 42 MU = 1,4
                CONE(MU,N)=CONE(MU,N)+AB*S(MU)
   42         CONTINUE
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
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO 43 MU = 1,12
                  CONE(MU,N)=CONE(MU,N)+AB*S(MU)
   43           CONTINUE
              END IF
            END IF
          END IF
        END IF
      END IF
   12 CONTINUE
      IF (ICOLBREI .EQ. 2) THEN
        DO 46 N = 1,NE2
          NU=NE1+2*(N-1)
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJA,JJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJA,JJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N.EQ.NE2) GO TO 46
          NUP1=NU+1
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJA,JJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJA,JJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJA,JJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJA,JJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
   46   CONTINUE
      END IF
      RETURN
   10 WRITE(99,100)KRA1
  100 FORMAT(5X,'ERRO IN EL53  PMGG RAGG KRA1=',I100)
      STOP
      END
