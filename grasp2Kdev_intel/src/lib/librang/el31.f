********************************************************************
*                                                                  *
      SUBROUTINE EL31(JJJA,JJJB,JA,JB,JJA,JJB,JJC,JJD,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 06  -------------  *
*                                                                  *
*     THIS PACKAGE EVALUATED THE CASES - 2111, 1211 ( + + - - ),   *
*     WHICH APPEARS IN CALCULATION MATRIX ELEMENTS BETWEEN         *
*     CONFIGURATIONS:                               N'1 = N1 - 1   *
*                                                   N'2 = N2 + 1   *
*                                                                  *
*     SUBROUTINE CALLED: COULOM,GG1222,ITREXG,IXJTIK,PERKO2,       *
*                        RECO,RECO2,SIXJ,SPEAK                     *
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
      DIMENSION J(2)
      DIMENSION S(12),IS(4),KAPS(4),KS(4)
      IF(NPEEL.LE.1)RETURN
      IIA=JLIST(JJA)
      IIB=JLIST(JJB)
      IIC=JLIST(JJC)
      IID=JLIST(JJD)
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
      CALL RECO2(JAA,JBB,J(2),0,IAT,REC)
      IF(IAT.EQ.0)RETURN
      IP1=ITREXG(J(1),J(1),J(1),J(2),IKK)+1
      IF(IKK.LE.0)RETURN
      IG1=IP1+IKK-1
      CALL RECO2(JAA,JBB,J(2),1,IAT,REC)
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
      END IF
C * * *                      * * *                      * * *
C     CASES 2111   + + - -        TRANSFORM TO  1112   + - - +
C           1211                                1112
C
      DO 2 I2=IP1,IG1,2
      KRA=(I2-1)/2
      IF (ICOLBREI .EQ. 1) THEN
        CALL COULOM(L2,L1,L1,L1,ID2(5),ID1(5),ID1(5),ID1(5),KRA,A1)
        IF(ABS(A1).LT.EPS)GO TO 2
      A1=-A1
      END IF
      AB=ZERO
      DO 3 I3=IP1,IG1,2
        J12=(I3-1)/2
        IFAZ=J(2)-J12+1
        IF((IFAZ/2)*2.NE.IFAZ)GO TO 3
        IF(IXJTIK(J(2),J(1),KRA*2,J(1),J(1),J12*2).EQ.0)GO TO 3
        CALL GG1222(IK2,IK1,BK2,BK1,ID2,ID1,BD2,
     *  BD1,J12,QM1,QM2,QM3,QM4,AA)
        IF(ABS(AA).LT.EPS) GO TO 3
        CALL SIXJ(J(2),J(1),KRA*2,J(1),J(1),J12*2,0,SI)
        AA=AA*SI*SQRT(DBLE(I3))
        IFAZ=2*J(1)+KRA*2+J12*2
        IF((IFAZ/4)*4.NE.IFAZ)AA=-AA
        AB=AB+AA
    3 CONTINUE
      AB=AB*REC
      IF(ABS(AB).LT.EPS)GO TO 2
C
C     TRANSFORM FANO & RACAH PHASE CONVENTION
C     TO CONDON & SHORTLEY PHASE CONVENTION
C
      IFAZFRCS = 1
      IFAZ=IK1(5)*IK1(4)+IK2(5)*IK2(4)-ID1(5)*ID1(4)-ID2(5)*ID2(4)
      IF((IFAZ/4)*4.NE.IFAZ)IFAZFRCS=-IFAZFRCS
C
      NN=0
      JB1=JBB-1
      DO 6 II=JAA,JB1
      IN=JLIST(II)
      NN=NQ1(IN)+NN
    6 CONTINUE
      IF((NN/2)*2.EQ.NN)AB=-AB
      IF (ICOLBREI .EQ. 1) THEN
         BB=A1*AB*DBLE(IFAZFRCS)
         CALL SPEAK(JJJA,JJJB,IIA,IIB,IIC,IID,KRA,BB)
      ELSE IF (ICOLBREI .EQ. 2) THEN
        N=(KRA-ND1)/2+1
        IF(((KRA-ND1)/2)*2 .EQ. (KRA-ND1)) THEN
          CALL CXK(S,IS,KAPS,KRA,KRA,2,1)
          IF(DABS(S(1)).GT.EPS) THEN
            BB=-S(1)*AB
            IF(DABS(BB).GT.EPS)
     :      CALL TALK(JJJA,JJJB,KRA,IS(1),IS(3),IS(2),IS(4),3,BB)
          END IF
        END IF
      END IF
    2 CONTINUE
      RETURN
      END
