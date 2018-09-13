************************************************************************
*                                                                      *
      SUBROUTINE START (IORB,ITYPE,P0,P,Q0,Q)
*                                                                      *
*   This subroutine sets up  P(1:6), Q(1:6),  required  to start the   *
*   integration for programs  OUT  and  SBSTEP .                       *
*                                                                      *
*   Arguments:                                                         *
*                                                                      *
*      IORB : (Input) Index of the orbital                             *
*      ITYPE: (Input) 1 = homogeneous equation; 2 = inhomogeneous      *
*             equation; 3 = variational equation                       *
*      P0   : (Input) Slope parameter                                  *
*      P    : (Output) P(1:6) are tabulated by this program            *
*      Q0   : (Output) First term in the series expansion of Q         *
*      Q    : (Output) Q(1:6) are tabulated by this program            *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTNXA,PNTNYA,PNTRXA,PNTRYA
      POINTER (PNTNXA,NXADUMMY)
      POINTER (PNTNYA,NYADUMMY)
      POINTER (PNTRXA,RXADUMMY)
      POINTER (PNTRYA,RYADUMMY)
      CHARACTER*2 NH
*
      POINTER (PNTRDA,DA(*))
      POINTER (PNTNDA,NDA(*))
*
      DIMENSION P(NNNP),Q(NNNP)
      DIMENSION RDP(6),RDQ(6),RSEP(6),RSEQ(6),SPEST(2:6),SQEST(2:6)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/CNC6/CNC6C(1:6,2:6)
     :      /DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF2/PNTRDA,PNTRXA,PNTRYA,
     :            PNTNDA,PNTNXA,PNTNYA,
     :            NDCOF,NXCOF,NYCOF,
     :            NDDIM,NXDIM,NYDIM
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      DATA MXITER/36/
*
*   Initialization
*
      OBC = 1.0D 00/C
      ZONC = Z*OBC
      GIORB = GAMA(IORB)
      KIORB = NAK(IORB)
      FKIORB = DBLE (KIORB)
      OMGI = 1.0D 00-GIORB
*
      OMGMK = OMGI-FKIORB
      OMGPK = OMGI+FKIORB
*
*   Determine P(1), Q(1): THESE STORE  R**(-GAMMA)*(P(1),Q(1));
*   set up  RSEP  and  RSEQ , the inhomogeneous terms
*
      IF ((ITYPE .EQ. 1) .OR. (ITYPE .EQ. 2)) THEN
         IF (NPARM .EQ. 0) THEN
            P(1) = P0
            IF (KIORB .LT. 0) THEN
               Q(1) = -P0*ZONC/(GIORB-FKIORB)
            ELSE
               Q(1) =  P0*(GIORB+FKIORB)/ZONC
            ENDIF
         ELSE
            IF (KIORB .LT. 0) THEN
               P(1) = P0
               Q(1) = 0.0D 00
            ELSE
               P(1) = 0.0D 00
               Q(1) = P0*(GIORB+FKIORB)/ZONC
            ENDIF
         ENDIF
         IF (ITYPE .EQ. 1) THEN
            DO 1 I = 1,6
               RSEP(I) = 0.0D 00
               RSEQ(I) = 0.0D 00
    1       CONTINUE
         ELSE
            RSEP1 = 0.0D 00
            RSEQ1 = 0.0D 00
            DO 2 I = 1,NDCOF
               JORB = NDA(I)
               PZERO = PZ(JORB)
               IF (NPARM .EQ. 0) THEN
                  P1 = PZERO
                  IF (KIORB .LT. 0) THEN
                     Q1 = -PZERO*ZONC/(GIORB-FKIORB)
                  ELSE
                     Q1 =  PZERO*(GIORB+FKIORB)/ZONC
                  ENDIF
                  SUMP =  ZONC*Q1
                  SUMQ = -ZONC*P1
               ELSE
                  IF (KIORB .LT. 0) THEN
                     P1 = PZERO
                     Q1 = 0.0D 00
                  ELSE
                     P1 = 0.0D 00
                     Q1 = PZERO*(GIORB+FKIORB)/ZONC
                  ENDIF
                  SUMP = 0.0D 00
                  SUMQ = 0.0D 00
               ENDIF
               FACTOR = DA(I)
               RSEP1 = RSEP1+FACTOR*(SUMP+OMGMK*P1)
               RSEQ1 = RSEQ1+FACTOR*(SUMQ+OMGPK*Q1)
    2       CONTINUE
            FACTOR = RP(1)
            RSEP(1) = FACTOR*RSEP1
            RSEQ(1) = FACTOR*RSEQ1
            DO 3 I = 2,6
               FACTOR = -RP(I)*R(I)**(-GIORB)
               RSEP(I) = FACTOR*XP(I)
               RSEQ(I) = FACTOR*XQ(I)
    3       CONTINUE
         ENDIF
      ELSEIF (ITYPE .EQ. 3) THEN
         P(1) = 0.0D 00
         Q(1) = 0.0D 00
         RSEP(1) = 0.0D 00
         RSEQ(1) = 0.0D 00
         DO 4 I = 2,6
            FACTOR = OBC*RP(I)*R(I)**(OMGI)
            RSEP(I) = -FACTOR*QF(I,IORB)
            RSEQ(I) =  FACTOR*PF(I,IORB)
    4    CONTINUE
      ENDIF
      Q0 = Q(1)
*
*   Set up  RDP  and  RDQ
*
      CSQ = C*C
      TWOCSQ = CSQ+CSQ
      ENERGY = E(IORB)
      ENEFAC = TWOCSQ-ENERGY
      DO 5 I = 1,6
         RI = R(I)
         RPI = RP(I)
         RIRPI = RI*RPI
         YPIRPI = YP(I)*RPI
         RDP(I) = -OBC*(ENEFAC*RIRPI+YPIRPI)
         RDQ(I) = -OBC*(ENERGY*RIRPI-YPIRPI)
    5 CONTINUE
*
*   Determine  P(2:6) , Q(2:6)
*
*   Initilizations for the iterations
*
      NITER = 0
      P1 = P(1)
      Q1 = Q(1)
      DIFMAW = MAX (ABS (P1),ABS (Q1))
*
      DO 6 I = 2,6
         P(I) = P1
         Q(I) = Q1
    6 CONTINUE
*
*   This is the largest factor by which any result will be
*   multiplied
*
      FACTOR = R(6)**GIORB
*
*   Now iterate
*
    7 NITER = NITER+1
      DIFMAX = 0.0D 00
      DO 9 J = 2,6
         SUMP = 0.0D 00
         SUMQ = 0.0D 00
         DO 8 I = 1,6
            COEFIJ = CNC6C(I,J)
            RPI = RP(I)
            PI = P(I)
            QI = Q(I)
            SUMP = SUMP+COEFIJ*(OMGMK*RPI*PI-RDP(I)*QI+RSEP(I))
            SUMQ = SUMQ+COEFIJ*(OMGPK*RPI*QI-RDQ(I)*PI+RSEQ(I))
    8    CONTINUE
         RJ = R(J)
         SUMP = SUMP/RJ
         SUMQ = SUMQ/RJ
         SPEST(J) = SUMP
         SQEST(J) = SUMQ
         DIFMAX = MAX (DIFMAX,ABS (SUMP-P(J)))
         DIFMAX = MAX (DIFMAX,ABS (SUMQ-Q(J)))
    9 CONTINUE
czou  IF (DIFMAX .LT. DIFMAW) THEN
         DO 10 I = 2,6
            P(I) = SPEST(I)
            Q(I) = SQEST(I)
   10    CONTINUE
         DIFMAW = DIFMAX
         DIFMAX = DIFMAX*FACTOR
         IF (DIFMAX .GT. ACCY) THEN
            IF (NITER .LT. MXITER) THEN
               GOTO 7
            ELSE
               WRITE (*,300) NP(IORB),NH(IORB),DIFMAX,NITER,ACCY
            ENDIF
         ENDIF
c     ELSE
c        DIFMAX = DIFMAX*FACTOR
c        IF (DIFMAX .GT. ACCY) THEN
c           WRITE (*,300) NP(IORB),NH(IORB),DIFMAX,NITER,ACCY
c        ENDIF
czou  ENDIF
c not convergent, using the initial P,Q
      IF (DIFMAX .GT. ACCY) THEN
        DO I = 2,6
          P(I) = P1
          Q(I) = Q1
        ENDDO   
      ENDIF
czou
*
*   All done
*
*   This is always true in GRASP2
*
      P(1) = 0.0D 00
      Q(1) = 0.0D 00
*
      DO 11 I = 2,6
         FACTOR = R(I)**GIORB
         P(I) = FACTOR*P(I)
         Q(I) = FACTOR*Q(I)
   11 CONTINUE
*
      RETURN
*
  300 FORMAT ('START: ',1I2,1A2,' subshell: accuracy ',1P,1D7.1,
     :       /' attained after ',1I2,' iterations; this fails the'
     :       /' accuracy criterion ',D7.1,'.')
*
      END
