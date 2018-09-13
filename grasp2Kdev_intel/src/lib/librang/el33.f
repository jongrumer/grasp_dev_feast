********************************************************************
*                                                                  *
      SUBROUTINE EL33(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD
     :                                                  ,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 08  -------------  *
*                                                                  *
*     THIS PACKAGE EVALUATED THE CASES - 2313, 3231, 3213, 2331    *
*                                                   ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 - 1   *
*                                                   N'2 = N2 + 1   *
*                                                                  *
*     SUBROUTINE CALLED: COULOM,EILE,GG1233,ITREXG,IXJTIK,         *
*                      JFAZE,PERKO2,RECO,REC3,SIXJ,SPEAK           *
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
      DIMENSION PMGG(30),RAGG(30),J(3)
      DIMENSION CONE(12,20),S(12),IS(4),KAPS(4),KS(4)
      IF(NPEEL.LE.1)RETURN
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT.EQ.0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      QM1=HALF
      QM2=-HALF
      QM3=HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      IF (ICOLBREI .EQ. 2) THEN
        IS(1)=IIA
        IS(2)=IIB
        IS(3)=IIC
        IS(4)=IID
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
      DO 4 I4=IP1,IG1,2
        KRA=(I4-1)/2
        KRA1=KRA+1
        IF(KRA1.GT.30)GO TO 10
        RAGG(KRA1)=ZERO
        PMGG(KRA1)=ZERO
        CALL REC3(JB,JA,JC,J(2),J(1),KRA*2,0,IAT,REC)
        IF(IAT.EQ.0)GO TO 4
        CALL GG1233(IK2,IK1,IK3,BK2,BK1,BK3,ID2,ID1,ID3,BD2,
     *  BD1,BD3,KRA,QM1,QM2,QM3,QM4,RAG)
        IF(ABS(RAG).LT.EPS) GO TO 4
        RAGG(KRA1)=RAG
        CALL REC3(JB,JA,JC,J(2),J(1),KRA*2,1,IAT,REC)
        PMGG(KRA1)=REC
    4 CONTINUE
      IFAZP=JFAZE(JB,JA,JC,JC)
C
C     TRANSFORM FANO & RACAH PHASE CONVENTION
C     TO CONDON & SHORTLEY PHASE CONVENTION
C
      IFAZFRCS=1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
      IF((IFAZ/4)*4.NE.IFAZ)IFAZFRCS=-IFAZFRCS
C
      IF(JA.GT.JB)GO TO 17
      JAA=JA
      JBB=JB
      GO TO 18
   17 JAA=JB
      JBB=JA
   18 NN=0
      JB1=JBB-1
      DO 16 II=JAA,JB1
      IN=JLIST(II)
      NN=NQ1(IN)+NN
   16 CONTINUE
      IF((NN/2)*2.EQ.NN)IFAZP=-IFAZP
      IF(IREZ.EQ.2)GO TO 5
C * * *                      * * *                      * * *
C     CASES 2313   + + - -        TRANSFORM TO  2133   + - + -
C           3231                                2133
C
    6 CONTINUE
      DO 1 I1=IP1,IG1,2
        KRA=(I1-1)/2
        KRA1=KRA+1
        IF(KRA1.GT.30)GO TO 10
        IF (ICOLBREI .EQ. 1) THEN
          CALL COULOM(L2,L3,L1,L3,ID2(5),ID3(5),ID1(5),ID3(5),KRA,A1)
          IF(ABS(A1).LT.EPS)GO TO 1
        END IF
        AA=PMGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 1
        AA=AA*RAGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 1
        AA=AA/SQRT(DBLE(I1))
        AA=AA*DBLE(IFAZP)
        IF (ICOLBREI .EQ. 1) THEN
          BB=A1*AA*DBLE(IFAZFRCS)
          CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
        ELSE IF (ICOLBREI .EQ. 2) THEN
          N=(KRA-ND1)/2+1
          CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
          IF(DABS(S(1)).GT.EPS) THEN
            BB=S(1)*AA
            IF(DABS(BB).GT.EPS)
     :    CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
          END IF
        END IF
    1 CONTINUE
      IF(IREZ.EQ.2)GO TO 7
C * * *                      * * *                      * * *
C     CASES 3213   + + - -        TRANSFORM TO  2133   + - + -
C           2331                                2133
C
    5 IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      DO 2 I2=IP2,IG2,2
      KRA=(I2-1)/2
      IF(KRA.GT.30)GO TO 10
      IF (ICOLBREI .EQ. 1) THEN
        CALL COULOM(L3,L2,L1,L3,ID3(5),ID2(5),ID1(5),ID3(5),KRA,A1)
        IF(ABS(A1).LT.EPS)GO TO 2
      END IF
      AB=ZERO
      DO 3 I3=IP1,IG1,2
        J12=(I3-1)/2
        KRA1=J12+1
        IF(KRA1.GT.30)GO TO 10
        AA=PMGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 3
        AA=AA*RAGG(KRA1)
        IF(ABS(AA).LT.EPS) GO TO 3
        IF(IXJTIK(J(1),J(3),KRA*2,J(3),J(2),J12*2).EQ.0)GO TO 3
        CALL SIXJ(J(1),J(3),KRA*2,J(3),J(2),J12*2,0,SI)
        AA=AA*SI*SQRT(DBLE(I3))
        AB=AB+AA
    3 CONTINUE
      IF(ABS(AB).LT.EPS)GO TO 2
      AB=AB*DBLE(IFAZP)
      IF (ICOLBREI .EQ. 1) THEN
        BB=A1*AB*DBLE(IFAZFRCS)
        CALL SPEAK(JJJA,JJJB,IIA,IIB,IID,IIC,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        NU=KRA
        IF(((NU-NE1)/2)*2 .EQ. (NU-NE1)) THEN
          IF((ITRIG(KS(1),KS(4),NU+NU+1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(3),NU+NU+1).NE.0)) THEN
            IF(NU .GT. 0) THEN
              N=(NU-NE1)/2+1
              CALL CXK(S,IS,KAPS,NU,KRA,1,2)
              DO 21 MU = 1,4
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
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO 22 MU = 1,4
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
                CALL CXK(S,IS,KAPS,NU,KRA,1,2)
                DO 23 MU = 1,12
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
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(4),IS(2),IS(3),1,CONE(1,N))
          CALL TALK(JJJA,JJJB,NU,IS(4),IS(1),IS(3),IS(2),1,CONE(2,N))
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(4),IS(3),IS(2),1,CONE(3,N))
          CALL TALK(JJJA,JJJB,NU,IS(4),IS(1),IS(2),IS(3),1,CONE(4,N))
          IF(N.EQ.NE2) GO TO 36
          NUP1=NU+1
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(4),IS(2),IS(3),2,CONE(5,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(3),IS(1),IS(4),2,CONE(6,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(1),IS(3),IS(2),2,CONE(7,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(2),IS(4),IS(1),2,CONE(8,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(4),IS(3),IS(2),2,CONE(9,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(2),IS(1),IS(4),2,CONE(10,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(1),IS(2),IS(3),2,CONE(11,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(3),IS(4),IS(1),2,CONE(12,N))
   36   CONTINUE
      END IF
      IF(IREZ.EQ.2)GO TO 6
    7 CONTINUE
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL33  PMGG RAGG')
      STOP
      END
