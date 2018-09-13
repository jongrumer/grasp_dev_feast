************************************************************************
*                                                                      *
      SUBROUTINE SOLVE (J,FAIL,INV,JP,NNP)
*                                                                      *
*   This subroutine  performs  step  2 in Algorithm 5.2 and 5.3 of C   *
*   Froese Fischer, Comput Phys Rep 3 (1986) 295. Some minor changes   *
*   have been made.                                                    *
*                                                                      *
*   Arguments:                                                         *
*                                                                      *
*      J     : (Input) The serial number of the orbital                *
*      JP    : (Output) The join point                                 *
*      FAIL  : (Output) If .true., the iterations did not yield an     *
*              acceptable solution (methods 1 and 2)                   *
*                                                                      *
*   Call(s) to: [RSCF92]: COUNT, ESTIM, IN, NEWE, OUT, cofpot,         *
*                         SETXUV, SETXV, SETXZ, START                  *
*               [LIB92]: DCBSRW, QUAD, SETPOT.                         *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 26 Sep 1993   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      CHARACTER*2 NH
      LOGICAL CHECK,FAIL,LDBPR,NOINVT
*
      DIMENSION PH(NNNP),QH(NNNP),PV(NNNP),QV(NNNP)
*
      POINTER (PNTRPF,PF(NNNP,1))
      POINTER (PNTRQF,QF(NNNP,1))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /INVT/NOINVT(NNNW)
     :      /NODE/NNODEP(NNNW)
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
     :      /SCF4/EPSMIN,EPSMAX,EMIN,EMAX,ZINF
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
*   Initialization
*
      CHECK = .NOT. NOINVT(J)
      FAIL = .FALSE.
      KOUNT = 0
      NAKJ = NAK(J)
cb alpha constant from lib/lib92/setcon.f
c     TWOCSQ = 2.0D 00*137.036**2
      TWOCSQ = 2.0D 00*C*C
*
*   Debug header
*
      IF (LDBPR(22)) WRITE (99,300) NP(J),NH(J)
*
      NLOOPS = MAX (NSOLV,3*NP(J))
*
      CALL ESTIM (J)
      ELAST = E(J)
      E(J) = EIGEN (J)
*
*   Checks on lower bounds
*
      IF (E(J) .LT. EPSMIN) THEN
         IF (E(J) .GT. 0.0D 00) THEN
            EPSMIN = E(J)
         ELSE
            WRITE (*,301) NP(J),NH(J),E(J)
            E(J) = EPSMIN
            IF ((ABS (EMIN-ELAST) .LE. 1.0D-06) .AND.
     :          (METHOD(J) .LE. 2)) THEN
               CALL DCBSRW (NP(J),NAKJ,ZINF,E(J),P0,P,Q,MTP0)
               WRITE (*,302) NP(J),NH(J),ZINF
               RETURN
            ENDIF
         ENDIF
      ENDIF
*
*   Check on upper bound
*
      IF (METHOD(J) .LE. 2) THEN
      IF (E(J) .GT. EPSMAX) THEN
         TENEJ = 10.0D 00*E(J)
         EPSMAX = MIN (TENEJ,TWOCSQ)
         EMAX = EPSMAX
         IF (E(J) .GT. TWOCSQ) THEN
            WRITE (*,303) NP(J),NH(J),E(J)
            E(J) = TWOCSQ
         ENDIF
      ENDIF
      ENDIF
*
*   Iteration loop begins here
*
    1 KOUNT = KOUNT+1
*
      P0 = PZ(J)
*
      IF (KOUNT .GT. 1) THEN
*
*   Check that bounds are ordered correctly
*
         IF (EPSMAX .LE. EPSMIN) WRITE (*,304) EPSMIN,EPSMAX,
     :                                           NP(J),NH(J),E(J)
*
         IF ((KOUNT .GT. NLOOPS) .OR.
     :       (EPSMAX-EPSMIN .LT.
     :          1.0D 00/(DBLE (NP(J))**3)))
     :      THEN
            WRITE (*,305) METHOD(J),NP(J),NH(J)
            WRITE (*,306) KOUNT-1,NLOOPS,
     :                      P0,
     :                      E(J),DELEPS,
     :                      EPSMIN,EPSMAX,
     :                      JP,MTP,
     :                      NNP,NNODEP(J),
     :                      SGN
            FAIL = .TRUE.
            RETURN
         ENDIF
      ENDIF
*
*   Set up arrays TF and TG; find join point
*
      CALL SETPOT (J,JP)
*
*   Set right-hand side to zero to form homogeneous equations;
*   integrate homogeneous equations outwards and inwards; store
*   small component at join point each time
*
      CALL SETXZ (J)
      ICASE = 1
      CALL START (J,ICASE,P0,PH,Q0,QH)
      CALL OUT (J,JP,PH,QH)
      QJPOH = QH(JP)
      CALL IN (J,JP,PH,QH,MTPH)
      QJPIH = QH(JP)
*
*   Set up right-hand side for inhomogeneous equations; integrate
*   inhomogeneous equations outwards and inwards; store small
*   component at join point each time
*
      CALL SETXUV (J)
      ICASE = 2
      CALL START (J,ICASE,P0,P,Q0,Q)
      CALL OUT (J,JP,P,Q)
      QJPOI = Q(JP)
      CALL IN (J,JP,P,Q,MTPI)
      QJPII = Q(JP)
*
*   Determine energy adjustment for methods 1 and 2
*
      IF (METHOD(J) .LE. 2) THEN
         TA(1) = 0.0D 00
         DO 2 I = 2,MTPI
            TA(I) = (P(I)**2+Q(I)**2)*RP(I)
    2    CONTINUE
         MTP = MTPI
         CALL QUAD (DNORM)
c        DELEPS = 137.036*P(JP)*(QJPII-QJPOI)/DNORM
         DELEPS = C*P(JP)*(QJPII-QJPOI)/DNORM
      ENDIF
*
*   Generate the continuous solution
*
      MTPC = MAX (MTPH,MTPI)
      ALFA = -(QJPII-QJPOI)/(QJPIH-QJPOH)
      P0H = P0
      P0 = P0*(1.0D 00+ALFA)
      DO 3 I = 1,MTPC
         P(I) = P(I)+ALFA*PH(I)
         Q(I) = Q(I)+ALFA*QH(I)
    3 CONTINUE
*
      IF ((METHOD(J) .EQ. 2) .OR. (METHOD(J) .EQ. 4)) THEN
*
*   Set up right-hand side for variational equations; integrate
*   variational equations outwards and inwards; store small
*   component at join point each time
*
         P0V = 0.0D 00
         CALL SETXV (J)
         ICASE = 3
         CALL START (J,ICASE,P0V,PV,Q0V,QV)
         CALL OUT (J,JP,PV,QV)
         QJPOV = QV(JP)
         CALL IN (J,JP,PV,QV,MTPV)
         QJPIV = QV(JP)
*
*   Generate continuous solutions
*
         MTPVC = MAX (MTPC,MTPV)
         ALFA = -(QJPIV-QJPOV)/(QJPIH-QJPOH)
         DO 4 I = 1,MTPVC
            PV(I) = PV(I)+ALFA*PH(I)
            QV(I) = QV(I)+ALFA*QH(I)
    4    CONTINUE
*
         TA(1) = 0.0D 00
         DO 5 I = 2,MTPC
            TA(I) = RP(I)*(P(I)**2+Q(I)**2)
    5    CONTINUE
         MTP = MTPC
         CALL QUAD (DNORM)
*
         MTP = MIN (MTPC,MTPVC)
         TA(1) = 0.0D 00
         DO 6 I = 2,MTP
            TA(I) = RP(I)*(P(I)*PV(I)+Q(I)*QV(I))
    6    CONTINUE
         CALL QUAD (CRNORM)
*
         TA(1) = 0.0D 00
         DO 7 I = 2,MTPVC
            TA(I) = RP(I)*(PV(I)**2+QV(I)**2)
    7    CONTINUE
         MTP = MTPVC
         CALL QUAD (DVNORM)
*
*   Determine deleps required to normalize new solution to
*   first order: modified form of solution to a quadratic
*   equation (see Press et al.)
*
         AA = DVNORM
         BB = CRNORM+CRNORM
         CC = DNORM-1.0D 00
         DISCR = BB*BB-4.0D 00*AA*CC
         IF (DISCR .GT. 0.0D 00) THEN
            QQ = -0.5D 00
     :           *(BB+SIGN (1.0D 00,BB)*SQRT (DISCR))
            ROOT1 = CC/QQ
            ROOT2 = QQ/AA
            PMX = 0.0D 00
            APM = 0.0D 00
            DO 8 I = 2,JP
               PATI = P(I)
               IF (PATI .GT. PMX) THEN
                  PMX = PATI
                  LOC1 = I
               ENDIF
               ABPI = ABS (PATI)
               IF (ABPI .GT. APM) THEN
                  APM = ABPI
                  LOC2 = I
               ENDIF
    8       CONTINUE
            IF (PMX .NE. 0.0D 00) THEN
               RATIO = APM/ABS (PMX)
               IF (RATIO .LT. 10.0D 00) THEN
                  LOC = LOC1
               ELSE
                  LOC = LOC2
               ENDIF
            ELSE
               LOC = LOC2
            ENDIF
            TEST1 = P(LOC)+ROOT1*PV(LOC)
            TEST2 = P(LOC)+ROOT2*PV(LOC)
            IF     ((TEST1 .GT. 0.0D 00) .AND.
     :              (TEST2 .LT. 0.0D 00)) THEN
               DELE = ROOT1
            ELSEIF ((TEST1 .LT. 0.0D 00) .AND.
     :              (TEST2 .GT. 0.0D 00)) THEN
               DELE = ROOT2
            ELSEIF ((TEST1 .GT. 0.0D 00) .AND.
     :              (TEST2 .GT. 0.0D 00)) THEN
               IF (TEST1 .LT. TEST2) THEN
                  DELE = ROOT1
               ELSE
                  DELE = ROOT2
               ENDIF
            ELSEIF ((TEST1 .LT. 0.0D 00) .AND.
     :              (TEST2 .LT. 0.0D 00)) THEN
               IF (TEST1 .GT. TEST2) THEN
                  DELE = ROOT1
               ELSE
                  DELE = ROOT2
               ENDIF
            ENDIF
         ELSE
            DELE = -BB/(AA+AA)
         ENDIF
*
*   Generate new solution
*
         MTP0 = MAX (MTPC,MTPVC)
         P0 = P0+DELE*ALFA*P0H
         DO 9 I = 2,MTP0
            P(I) = P(I)+DELE*PV(I)
            Q(I) = Q(I)+DELE*QV(I)
    9    CONTINUE
      ELSE
         MTP0 = MTPC
      ENDIF
*
*   Debug printout
*
      IF (LDBPR(23)) CALL PRWF (J)
*
*   Count nodes in large component function; determine sign at first
*   oscillation, effective quantum number; note that node counting
*   is never enforced on the small component
*
      CALL COUNT (P,MTP0,NNP,SGN)
*
*   DEBUG PRINTOUT
*
      IF (LDBPR(22))
     :   WRITE (99,306) KOUNT,NLOOPS,
     :                   P0,
     :                   E(J),DELEPS,
     :                   EPSMIN,EPSMAX,
     :                   JP,MTP,
     :                   NNP,NNODEP(J),
     :                   SGN
*
*   Proceed according to method
*
      IF (METHOD(J) .GT. 2) THEN
         IF (CHECK .AND. (SGN .LT. 0.0D 00)) THEN
            INV = 1
            P0 = -P0
            DO 10 I = 2,MTP0
               P(I) = -P(I)
               Q(I) = -Q(I)
   10       CONTINUE
         ENDIF
      ELSE
         MX = NNP-NNODEP(J)
         NPRIME = NNP+NKL(J)+1
         CALL NEWE (J,SGN,NPRIME,MX,DELEPS,FAIL,INV)
         IF (FAIL) GOTO 1
      ENDIF
*
*   Solution found
*
      RETURN
*
  300 FORMAT (//' Debug printout active; orbital: ',1I2,1A2)
  301 FORMAT (' E(',1I2,1A2,') = ',1PD11.4,'; adjusted to EPSMIN')
  302 FORMAT (' Returned hydrogenic function for ',1I2,1A2,' with',
     :        ' effective charge ',F7.3)
  303 FORMAT (' E(',1I2,1A2,') = ',1PD11.4,'; adjusted to TWOCSQ')
  304 FORMAT (' Warning: difficulty with node-counting procedure'
     :       /' lower bound on energy (',1PD11.4,') exceeds upper'
     :       ,' bound (',1D11.4,'; E(',1I2,1A2,') = ',1D11.4)
  305 FORMAT (' Method ',1I1,' unable to solve for ',1I2,1A2,' orbital')
  306 FORMAT (' Iteration number: ',1I2,', limit: ',1I2
     :       /' Present estimate of P0; ',1D21.14
     :       /' Present estimate of E(J): ',1D21.14,', DELEPS: ',1D21.14
     :       /' Lower bound on energy: ',1D21.14,', upper bound: '
     :       ,1D21.14
     :       /' Join point: ',1I4,', Maximum tabulation point:',1I4
     :       /' Number of nodes counted: ',1I2,', Correct number: ',1I2
     :       /' Sign of P at first oscillation: ',F3.0)
*
      END
