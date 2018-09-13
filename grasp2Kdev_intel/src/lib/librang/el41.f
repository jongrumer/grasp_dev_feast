********************************************************************
*                                                                  *
      SUBROUTINE EL41(JJJA,JJJB,JA,JB,JC,IREZ,JJA,JJB,JJC,JJD,
     :                                                    ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 10  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 + 1        *
*                                              N'2 = N2 + 1        *
*                                              N'3 = N3 - 2,       *
*     WHEN IREZ = 1   . . . . . . . . . . . . . . . . . . .        *
*                                              N'1 = N1 - 1        *
*                                              N'2 = N2 - 1        *
*                                              N'3 = N3 + 2,       *
*     WHEN IREZ = 2   . . . . . . . . . . . . . . . . . . .        *
*                                                                  *
*     SUBROUTINE CALLED: COULOM,EILE,GG1233,ITREXG,IXJTIK,         *
*                        JFAZE,PERKO2,RECO,REC3,SIXJ,SPEAK         *
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
      DIMENSION J(3)
      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
      CALL EILE(JA,JB,JC,JAA,JBB,JCC)
      IF(NPEEL.LE.1)RETURN
      CALL RECO(JAA,JCC,JBB,JBB,2,IAT)
      IF(IAT.EQ.0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      IC=JLIST(JC)
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
      IF(IREZ.EQ.1)GO TO 20
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      GO TO 21
   20 QM1=-HALF
      QM2=-HALF
      QM3=HALF
      QM4=HALF
   21 CALL PERKO2(JA,JB,JC,JA,3)
      J(1)=ID1(3)
      J(2)=ID2(3)
      J(3)=ID3(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      L3=(J(3)+1)/2
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
        IF(IBRD .LE. 0)RETURN
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
    8   CONTINUE
      END IF
      IFAZP=JFAZE(JC,JA,JB,JC)
      IFAZFRCS = 1
C
C     TRANSFORM FANO & RACAH PHASE CONVENTION
C     TO CONDON & SHORTLEY PHASE CONVENTION
C
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)+
     *IK3(5)*IK3(4)-ID3(5)*ID3(4)
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
C * * *                      * * *                      * * *
C     CASES 3312   + + - -        TRANSFORM TO  1233   - - + +
C           3321                                1233
C                                                    (IREZ = 1)
C     OR
C     CASES 1233   + + - -        TRANSFORM TO  1233   + + - -
C           2133                                1233
C                                                    (IREZ = 2)
      IP1=ITREXG(J(2),J(1),J(3),J(3),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      IP2=ITREXG(J(3),J(1),J(2),J(3),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      DO 2 I2=IP2,IG2,2
      KRA=(I2-1)/2
C
      IF (ICOLBREI .EQ. 1) THEN
        IF(IREZ.EQ.2)GO TO 19
        CALL COULOM(L3,L3,L1,L2,ID3(5),ID3(5),ID1(5),ID2(5),KRA,A1)
        GO TO 22
   19   CALL COULOM(L1,L2,L3,L3,ID1(5),ID2(5),ID3(5),ID3(5),KRA,A1)
   22   IF(ABS(A1).LT.EPS)GO TO 2
      END IF
C
      AB=ZERO
      DO 3 I3=IP1,IG1,2
      J12=(I3-1)/2
      IFAZ=J(2)-J12+1
      IF(IREZ.EQ.2)IFAZ=J(1)-J12+1
      IF((IFAZ/2)*2.NE.IFAZ)GO TO 3
      CALL REC3(JA,JB,JC,J(1),J(2),J12*2,0,IAT,AA)
      IF(IAT.EQ.0)GO TO 3
      IF(IXJTIK(J(3),J(1),KRA*2,J(2),J(3),J12*2).EQ.0)GO TO 3
      CALL GG1233(IK1,IK2,IK3,BK1,BK2,BK3,ID1,ID2,ID3,BD1,
     *BD2,BD3,J12,QM1,QM2,QM3,QM4,AA)
      IF(ABS(AA).LT.EPS) GO TO 3
      CALL REC3(JA,JB,JC,J(1),J(2),J12*2,1,IAT,REC)
      AA=AA*REC
      CALL SIXJ(J(3),J(1),KRA*2,J(2),J(3),J12*2,0,SI)
      AA=AA*SI*SQRT(DBLE(I3))
      IFAZ=J(3)+J(1)+2*J12+2*KRA
      IF(IREZ.EQ.2)IFAZ=J(2)+J(3)+2*J12+2*KRA
      IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
      AB=AB+AA
    3 CONTINUE
      IF(ABS(AB).LT.EPS)GO TO 2
      AB=-AB*DBLE(IFAZP)
      IF (ICOLBREI .EQ. 1) THEN
        BB=A1*AB*DBLE(IFAZFRCS)
        CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
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
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(2),IS(4),1,COND(1,N))
          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(4),IS(2),1,COND(2,N))
          CALL TALK(JJJA,JJJB,NU,IS(1),IS(3),IS(4),IS(2),1,COND(3,N))
          CALL TALK(JJJA,JJJB,NU,IS(3),IS(1),IS(2),IS(4),1,COND(4,N))
          IF(N.EQ.ND2) GO TO 36
          NUP1=NU+1
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(2),IS(4),2,COND(5,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(1),IS(3),2,COND(6,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(4),IS(2),2,COND(7,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(3),IS(1),2,COND(8,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(1),IS(3),IS(4),IS(2),2,COND(9,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(4),IS(2),IS(1),IS(3),2,COND(10,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(3),IS(1),IS(2),IS(4),2,COND(11,N))
          CALL TALK(JJJA,JJJB,NUP1,IS(2),IS(4),IS(3),IS(1),2,COND(12,N))
   36   CONTINUE
      END IF
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL41  PMGG RAGG')
      STOP
      END
