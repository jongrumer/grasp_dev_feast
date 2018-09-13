********************************************************************
*                                                                  *
      SUBROUTINE EL2(JJA,JJB,JA,JB,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 04  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 2        *
*                                              N'2 = N2 + 2        *
*                                                                  *
*     SUBROUTINE CALLED: COULOM,GG1122,ITREXG,IXJTIK,PERKO2,       *
*                        RECO,RECO2,SIXJ,SPEAK                     *
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
      DIMENSION J(2)
      DIMENSION COND(12,20),S(12),IS(4),KAPS(4),KS(4)
      IF(NPEEL.LE.1)RETURN
      IF(JA.GT.JB)GO TO 4
      JAA=JA
      JBB=JB
      GO TO 1
    4 JAA=JB
      JBB=JA
    1 CALL RECO(JAA,JBB,JBB,JBB,1,IAT)
      IF(IAT.EQ.0)RETURN
      IA=JLIST(JA)
      IB=JLIST(JB)
      QM1=HALF
      QM2=HALF
      QM3=-HALF
      QM4=-HALF
      CALL PERKO2(JA,JB,JA,JA,2)
      J(1)=ID1(3)
      J(2)=ID2(3)
      L1=(J(1)+1)/2
      L2=(J(2)+1)/2
      IP1=ITREXG(J(1),J(1),J(2),J(2),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
C * * *                      * * *                      * * *
C     THE CASE 1122   + + - -
C
      IP2=ITREXG(J(1),J(2),J(1),J(2),IKK)+1
      IF(IKK.LE.0) RETURN
      IG2=IP2+IKK-1
      IF (ICOLBREI .EQ. 2) THEN
        IS(1)=IA
        IS(2)=IA
        IS(3)=IB
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
      DO 2 I2=IP2,IG2,2
      KRA=(I2-1)/2
      IF (ICOLBREI .EQ. 1) THEN
        CALL COULOM(L1,L1,L2,L2,ID1(5),ID1(5),ID2(5),ID2(5),KRA,A1)
        IF(ABS(A1).LT.EPS)GO TO 2
        A1=-HALF*A1
      END IF
      AB=ZERO
      DO 3 I3=IP1,IG1,2
        J12=(I3-1)/2
        CALL RECO2(JAA,JBB,J12*2,0,IAT,REC)
        IF(IAT.NE.0) THEN
          IF(IXJTIK(J(1),J(2),KRA*2,J(2),J(1),J12*2).NE.0) THEN
            CALL GG1122(J12,J12,QM1,QM2,QM3,QM4,AA)
            IF(ABS(AA).GT.EPS) THEN
              CALL RECO2(JAA,JBB,J12*2,1,IAT,REC)
              AA=AA*REC
              CALL SIXJ(J(1),J(2),KRA*2,J(2),J(1),J12*2,0,SI)
              AA=AA*SI*SQRT(DBLE(I3))
              IFAZ=IK1(3)+IK2(3)+KRA*2+J12*2
              IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
              AB=AB+AA
           END IF
         END IF
       END IF
    3 CONTINUE
C
C     TRANSFORM FANO & RACAH PHASE CONVENTION
C     TO CONDON & SHORTLEY PHASE CONVENTION
C
      IFAZFRCS=1
      IFAZ1=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
      IF((IFAZ1/4)*4.NE.IFAZ1)IFAZFRCS=-IFAZFRCS
C
      IF (ICOLBREI .EQ. 1) THEN
        BB=A1*AB*DBLE(IFAZFRCS)
        IF(ABS(BB).GT.EPS)CALL SPEAK(JJA,JJB,IA,IA,IB,IB,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        NU=KRA
        IF(((NU-ND1)/2)*2 .EQ. (NU-ND1)) THEN
          IF((ITRIG(KS(1),KS(3),NU+NU+1).NE.0) .AND.
     :       (ITRIG(KS(2),KS(4),NU+NU+1).NE.0)) THEN
            N=(NU-ND1)/2+1
            IF(NU .GT. 0) THEN
              CALL CXK(S,IS,KAPS,NU,KRA,1,1)
              DO 21 MU = 1,4
                COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
   21         CONTINUE
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
              DO 22 MU = 1,4
                COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
   22         CONTINUE
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
                DO 23 MU = 1,12
                  COND(MU,N)=COND(MU,N)-HALF*AB*S(MU)
   23           CONTINUE
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
      RETURN
      END
