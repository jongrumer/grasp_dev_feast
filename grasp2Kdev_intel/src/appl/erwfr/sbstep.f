************************************************************************
*                                                                      *
      SUBROUTINE SBSTEP (IORB,NSTRT,NEND,P,Q)
*                                                                      *
*   This  subroutine continues the solution of the homogeneous Dirac   *
*   radial equation  from tabulation point NSTRT to tabulation point   *
*   NEND. The algorithm of J E Sienkiewicz and W E Baylis, J Phys B:   *
*   At Mol Phys 20 (1987) 5145-5156, p 5155, is used.                  *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL OUT
*
      DIMENSION P(NNNP),Q(NNNP)
*
      COMMON/DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /SBC/C1,C2,C3,C4,C5,C6
*
      PARAMETER (B1 = -0.50D 00,
     :           B2 =  0.50D 00,
     :           B3 =  0.75D 00,
     :           B4 =  0.25D 00)
*
*   Initializations
*
      TBH = 2.0D 00/H
      TC1 = 2.0D 00*C1
*
*   Determine whether integration is inward or outwarD
*
      IDIFF = NEND-NSTRT
      IF (IDIFF .GT. 0) THEN
         OUT = .TRUE.
      ELSEIF (IDIFF .EQ. 0) THEN
         RETURN
      ELSEIF (IDIFF .LT. 0) THEN
         OUT = .FALSE.
      ENDIF
*
*   Overall initializations
*
      FK = DBLE (NAK(IORB))
      CC = C1*FK*H
*
*   Perform integration depending on case
*
      IF (OUT) THEN
*
*   Initialization for outward integration
*
         LOC = NSTRT
         PJ = P(LOC)
         QJ = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJ = -FAC*PJ-TBH*TF(LOC)*QJ
         QPJ =  FAC*QJ-TBH*TG(LOC)*PJ
*
         LOC = LOC-1
         PJM1 = P(LOC)
         QJM1 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM1 = -FAC*PJM1-TBH*TF(LOC)*QJM1
         QPJM1 =  FAC*QJM1-TBH*TG(LOC)*PJM1
*
         LOC = LOC-1
         PJM2 = P(LOC)
         QJM2 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM2 = -FAC*PJM2-TBH*TF(LOC)*QJM2
         QPJM2 =  FAC*QJM2-TBH*TG(LOC)*PJM2
*
         LOC = LOC-1
         PJM3 = P(LOC)
         QJM3 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM3 = -FAC*PJM3-TBH*TF(LOC)*QJM3
         QPJM3 =  FAC*QJM3-TBH*TG(LOC)*PJM3
*
         LOC = LOC-1
         PJM4 = P(LOC)
         QJM4 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJM4 = -FAC*PJM4-TBH*TF(LOC)*QJM4
         QPJM4 =  FAC*QJM4-TBH*TG(LOC)*PJM4
*
*   March out
*
         J = NSTRT-1
    1    J = J+1
*
            RPPJ =   B1*PJ  + B2*PJM1 + B3*PJM2 + B4*PJM3
     :             + C2*PPJ  + C3*PPJM1 + C4*PPJM2 + C5*PPJM3 + C6*PPJM4
            RQPJ =   B1*QJ  + B2*QJM1 + B3*QJM2 + B4*QJM3
     :             + C2*QPJ  + C3*QPJM1 + C4*QPJM2 + C5*QPJM3 + C6*QPJM4
*
            JP1 = J+1
            CCRPOR = CC*RPOR(JP1)
            CPJP1 = 1.0D 00+CCRPOR
            CMJP1 = 1.0D 00-CCRPOR
            FJP1 = TC1*TF(JP1)
            GJP1 = TC1*TG(JP1)
            DENOM = CPJP1*CMJP1-GJP1*FJP1
            FACTOR = 1.0D 00/DENOM
            PJP1 = (CMJP1*RPPJ-FJP1*RQPJ)*FACTOR
            QJP1 = (CPJP1*RQPJ-GJP1*RPPJ)*FACTOR
            P(JP1) = PJP1
            Q(JP1) = QJP1
*
            IF (JP1 .LT. NEND) THEN
*
               PPJM4 = PPJM3
               QPJM4 = QPJM3
*
               PJM3 = PJM2
               QJM3 = QJM2
               PPJM3 = PPJM2
               QPJM3 = QPJM2
*
               PJM2 = PJM1
               QJM2 = QJM1
               PPJM2 = PPJM1
               QPJM2 = QPJM1
*
               PJM1 = PJ
               QJM1 = QJ
               PPJM1 = PPJ
               QPJM1 = QPJ
*
               PJ = PJP1
               QJ = QJP1
               FAC = FK*RPOR(JP1)
               PPJ = -FAC*PJ-TBH*TF(JP1)*QJ
               QPJ =  FAC*QJ-TBH*TG(JP1)*PJ
*
               GOTO 1
*
            ENDIF
      ELSE
*
*   Initializations for inward integration
*
         LOC = NSTRT
         PJ = P(LOC)
         QJ = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJ = -FAC*PJ-TBH*TF(LOC)*QJ
         QPJ =  FAC*QJ-TBH*TG(LOC)*PJ
*
         LOC = LOC+1
         PJP1 = P(LOC)
         QJP1 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP1 = -FAC*PJP1-TBH*TF(LOC)*QJP1
         QPJP1 =  FAC*QJP1-TBH*TG(LOC)*PJP1
*
         LOC = LOC+1
         PJP2 = P(LOC)
         QJP2 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP2 = -FAC*PJP2-TBH*TF(LOC)*QJP2
         QPJP2 =  FAC*QJP2-TBH*TG(LOC)*PJP2
*
         LOC = LOC+1
         PJP3 = P(LOC)
         QJP3 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP3 = -FAC*PJP3-TBH*TF(LOC)*QJP3
         QPJP3 =  FAC*QJP3-TBH*TG(LOC)*PJP3
*
         LOC = LOC+1
         PJP4 = P(LOC)
         QJP4 = Q(LOC)
         FAC = FK*RPOR(LOC)
         PPJP4 = -FAC*PJP4-TBH*TF(LOC)*QJP4
         QPJP4 =  FAC*QJP4-TBH*TG(LOC)*PJP4
*
*   March in
*
         J = NSTRT+1
    2    J = J-1
*
            RPMJ =   B1*PJ  + B2*PJP1 + B3*PJP2 + B4*PJP3
     :             - C2*PPJ  - C3*PPJP1 - C4*PPJP2 - C5*PPJP3 - C6*PPJP4
            RQMJ =   B1*QJ  + B2*QJP1 + B3*QJP2 + B4*QJP3
     :             - C2*QPJ  - C3*QPJP1 - C4*QPJP2 - C5*QPJP3 - C6*QPJP4
*
            JM1 = J-1
            CCRPOR = CC*RPOR(JM1)
            CPJM1 = 1.0D 00+CCRPOR
            CMJM1 = 1.0D 00-CCRPOR
            FJM1 = TC1*TF(JM1)
            GJM1 = TC1*TG(JM1)
            DENOM = CPJM1*CMJM1-GJM1*FJM1
            FACTOR = 1.0D 00/DENOM
            PJM1 = (CPJM1*RPMJ+FJM1*RQMJ)*FACTOR
            QJM1 = (CMJM1*RQMJ+GJM1*RPMJ)*FACTOR
            P(JM1) = PJM1
            Q(JM1) = QJM1
*
            IF (JM1 .GT. NEND) THEN
*
               PPJP4 = PPJP3
               QPJP4 = QPJP3
*
               PJP3 = PJP2
               QJP3 = QJP2
               PPJP3 = PPJP2
               QPJP3 = QPJP2
*
               PJP2 = PJP1
               QJP2 = QJP1
               PPJP2 = PPJP1
               QPJP2 = QPJP1
*
               PJP1 = PJ
               QJP1 = QJ
               PPJP1 = PPJ
               QPJP1 = QPJ
*
               PJ = PJM1
               QJ = QJM1
               FAC = FK*RPOR(JM1)
               PPJ = -FAC*PJ-TBH*TF(JM1)*QJ
               QPJ =  FAC*QJ-TBH*TG(JM1)*PJ
*
               GOTO 2
*
            ENDIF
      ENDIF
*
      RETURN
      END
