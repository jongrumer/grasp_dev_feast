************************************************************************
*                                                                      *
      SUBROUTINE SOLVH (IORB,FAIL)
*                                                                      *
*   This routine solves the homogeneous Dirac radial equation.         *
*                                                                      *
*   Arguments:  IORB : (Input) Index of orbital                        *
*               FAIL : (Output) .TRUE. if solution not obtained        *
*                                                                      *
*   The direct potential is assumed tabulated in the COMMON array YP   *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*               [RSCF92]: COUNT, SBSTEP, SETPOT, START, TAIL.          *
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
      LOGICAL FAIL
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      DATA MXK /75/
*
*   A solution is deemed continuous when the relative mismatch in the
*   small component is within EPSLON
*
      EPSLON = ACCY*0.1E 00
*
*   Establish the number of nodes in the large component
*
      NPIORB = NP(IORB)
      NKIORB = NAK(IORB)
      NAKABS = ABS (NKIORB)
      IF (NKIORB .LT. 0) THEN
         NLIORB = NAKABS-1
      ELSE
         NLIORB = NAKABS
      ENDIF
      NNP = NPIORB-NLIORB-1
*
*   Establish the bounds on, and an estimate of, the eigenvalue
*
      CSQ = C*C
      ALPHA = 1.0D 00/C
      NREL = NPIORB-NAKABS
      FKABS = DBLE (NAKABS)
      FKAP2 = FKABS*FKABS
*
      ZALPHA = YP(N)*ALPHA
      IF (ZALPHA .LT. FKABS) THEN
         GAMMA = SQRT (FKAP2-ZALPHA*ZALPHA)
         EBYM = 1.0D 00
     :          /SQRT ( 1.0D 00
     :                 +(ZALPHA/(GAMMA+NREL+0.5D 00))**2)
         EMIN = (1.0D00-EBYM)*CSQ
      ELSE
         EMIN = 0.25D 00*CSQ/DBLE (NPIORB*NPIORB)
      ENDIF
*
      ZALPHA = Z*ALPHA
*
      IF (ZALPHA .LT. FKABS) THEN
         GAMMA = SQRT (FKAP2-ZALPHA*ZALPHA)
         EBYM = 1.0D 00
     :          /SQRT ( 1.0D 00
     :                 +(ZALPHA/(GAMMA+NREL))**2)
         E(IORB) = (1.0D 00-EBYM)*CSQ
      ELSE
         E(IORB) = CSQ
      ENDIF
*
      IF (ZALPHA .LT. FKABS) THEN
         GAMMA = SQRT (FKAP2-ZALPHA*ZALPHA)
         EBYM = 1.0D 00
     :          /SQRT ( 1.0D 00
     :                 +(ZALPHA/(GAMMA+NREL-0.5D 00))**2)
         EMAX = (1.0D 00-EBYM)*CSQ
      ELSE
         EMAX = CSQ+CSQ
      ENDIF
*
      DELE = 0.0D 00
*
*   Initialize
*
      FAIL = .FALSE.
      KOUNT = -1
*
*   Iteration loop begins here
*
    1 KOUNT = KOUNT+1
      IF (KOUNT .GT. MXK) THEN
         FAIL = .TRUE.
         RETURN
      ENDIF
*
*   Generate estimate of eigenvalue for this iteration
*
      EEST = E(IORB)+DELE
      IF ((EEST .GT. EMIN) .AND. (EEST .LT. EMAX)) THEN
         E(IORB) = EEST
      ELSE
         E(IORB) = 0.5D 00*(EMIN+EMAX)
      ENDIF
*
*   Set up arrays TF and TG; find join point
*
      CALL SETPOT (IORB,JP)
*
*   Initialize outward integration
*
      ITYPE = 1
      CALL START (IORB,ITYPE,PZ(IORB),PF(1,IORB),Q0,QF(1,IORB))
*
*   Continue outward integration to classical turning point
*
      NSTRT = 6
      CALL SBSTEP (IORB,NSTRT,JP,PF(1,IORB),QF(1,IORB))
      PFJPO = PF(JP,IORB)
      QFJPO = QF(JP,IORB)
*
*   Initialize inward integration
*
      CALL TAIL (IORB,PF(1,IORB),QF(1,IORB),JP,MTP)
*
*   Continue inward integration to classical turning point
*
      NSTRT = MTP
      CALL SBSTEP (IORB,NSTRT,JP,PF(1,IORB),QF(1,IORB))
      PFJPI = PF(JP,IORB)
      QFJPI = QF(JP,IORB)
*
*   Make large component continuous, determine mismatch in small
*   component
*
      RATIO = PFJPO/PFJPI
      DO 2 I = JP,N
         PF(I,IORB) = PF(I,IORB)*RATIO
         QF(I,IORB) = QF(I,IORB)*RATIO
    2 CONTINUE
*
*   Count nodes
*
      CALL COUNT (PF(1,IORB),MTP,NPC,SGN)
*
*   Correct the energy estimate if the number of nodes is wrong
*
      IF (NPC .GT. NNP) THEN
         EMIN =  E(IORB)
         DELE =  0.5D 00*(EMAX-EMIN)
         GOTO 1
      ELSEIF (NPC .LT. NNP) THEN
         EMAX =  E(IORB)
         DELE = -0.5D 00*(EMAX-EMIN)
         GOTO 1
      ENDIF
*
*   Correct number of nodes
*
*   Compute 'norm' of solution
*
      TA(1) = 0.0D 00
      DO 3 I = 2,N
         TA(I) = (PF(I,IORB)**2+QF(I,IORB)**2)*RP(I)
    3 CONTINUE
      MTP = N
      CALL QUAD (DNORM)
*
*   Determine correction to eigenvalue from magic formula
*   correct slope at origin
*
      QFJPI = QFJPI*RATIO
      DMSMCH = QFJPI-QFJPO
      IF (ABS (DMSMCH/QFJPO) .GT. EPSLON) THEN
         DELE = C*PF(JP,IORB)*DMSMCH/DNORM
         IF (DELE .LT. 0.0D 00) THEN
            EMAX = EMAX*(1.0D 00
     :                  -0.2E 00*ABS (DELE/E(IORB)))
         ELSE
            EMIN = EMIN*(1.0D 00
     :                  +0.2E 00*ABS (DELE/E(IORB)))
         ENDIF
         PZ(IORB) = PZ(IORB)/SQRT (DNORM)
         GOTO 1
      ENDIF
*
*   Normalize
*
      DNFAC = 1.0D 00/SQRT (DNORM)
      PZ(IORB) = PZ(IORB)*DNFAC
      DO 4 I = 1,N
         PF(I,IORB) = PF(I,IORB)*DNFAC
         QF(I,IORB) = QF(I,IORB)*DNFAC
    4 CONTINUE
*
*   Find maximum tabulation point
*
      I = N+1
    5 I = I-1
      IF (ABS (PF(I,IORB)) .LT. EPSLON) THEN
         PF(I,IORB) = 0.0D 00
         QF(I,IORB) = 0.0D 00
         GOTO 5
      ELSE
         MF(IORB) = I
      ENDIF
*
      RETURN
      END
