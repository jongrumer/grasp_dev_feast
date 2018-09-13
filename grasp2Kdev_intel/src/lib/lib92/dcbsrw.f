************************************************************************
*                                                                      *
      SUBROUTINE DCBSRW (N,KAPPA,Z,E,RG0,RG,RF,MTP)
*                                                                      *
*   This subroutine computes the  Dirac-Coulomb  bound-state orbital   *
*   radial wavefunction.   Equations (13.5) and (13.5') of  Akhiezer   *
*   and Berestetskii modified to ensure positive slope at the origin   *
*   for RG are used.                                                   *
*                                                                      *
*   The arguments are as follows:                                      *
*                                                                      *
*      N    : (Input) The (usual) principal quantum number             *
*      KAPPA: (Input) The relativistic angular quantum number          *
*      Z    : (Input) The effective nuclear charge                     *
*      E    : (Output) The Dirac-Coulomb Eigenenergy                   *
*      RG0  : (Output) Coefficient of the leading term in the          *
*                      series expansion of the large component         *
*                      near the origin                                 *
*      RG   : (Output) r times the large component wavefunction of     *
*                      Akhiezer and Berestetskii                       *
*      RF   : (Output) r times the small component wavefunction of     *
*                      Akhiezer and Berestetskii                       *
*      MTP  : (Output) Maximum tabulation point                        *
*                                                                      *
*   Call(s) to: [LIB92]: CGAMMA.                                       *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last Update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      DIMENSION RG(NNNP),RF(NNNP)
*
      COMMON/DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,NTP
     :      /TATB/TA(NNN1),TB(NNN1),MTPAB
*
*   Ensure that the principal quantum number is physical
*
      IF (N .LE. 0) THEN
         WRITE (*,300)
         WRITE (*,301) N
         STOP
      ENDIF
*
*   Ensure that the angular quantum number is physical
*
      IF     (     KAPPA  .EQ. 0) THEN
         WRITE (*,300)
         WRITE (*,302)
         STOP
      ELSEIF (     KAPPA  .EQ. N) THEN
         WRITE (*,300)
         WRITE (*,303) KAPPA,N
         STOP
      ELSEIF (ABS (KAPPA) .GT. N) THEN
         WRITE (*,300)
         WRITE (*,303) KAPPA,N
         STOP
      ENDIF
*
*   Ensure that the charge is physical
*
      IF     (Z .LE. 0.0D 00) THEN
         WRITE (*,300)
         WRITE (*,304) Z
         STOP
      ELSEIF (Z .GT. C) THEN
         WRITE (*,300)
         WRITE (*,305) Z,C
         STOP
      ENDIF
*
*   Atomic units
*
      ALFA = 1.0D 00/C
*
*   Now determine all the parameters
*
      FN = DBLE (N)
      FKAPPA = DBLE (KAPPA)
      K = ABS (KAPPA)
      FK = DBLE (K)
      NR = N-K
      FNR = DBLE (NR)
      ZALFA = Z*ALFA
      GAMMA = SQRT (FK*FK-ZALFA*ZALFA)
      TWOGP1 = GAMMA+GAMMA+1.0D 00
      BIGN = SQRT (FN*FN-2.0D 00*FNR*(FK-GAMMA))
      EPS = 1.0D 00
     :      /SQRT (1.0D 00+(ZALFA/(GAMMA+FNR))**2)
*
*   EPS is the total energy divided by C*C; this must be converted
*   to the units and reference energy of GRASP
*
      E = (1.0D 00-EPS)*C*C
*
*   Now the normalization constants
*
      NRFAC = 1
      DO 1 I = 1,NR
         NRFAC = NRFAC*I
    1 CONTINUE
*
      ARGR = TWOGP1+FNR
      ARGI = 0.0D 00
      CALL CGAMMA (ARGR,ARGI,RGAMM1,DUMMY)
      ARGR = TWOGP1
      CALL CGAMMA (ARGR,ARGI,RGAMM2,DUMMY)
*
      FAC = - SQRT (RGAMM1)/
     :     (RGAMM2*SQRT (DBLE (NRFAC)))*
     :      SQRT (Z/(2.0D 00*BIGN*BIGN*(BIGN-FKAPPA)))
*
*   Ensure that the slope of the large-component function is
*   positive at the origin
*
      IF (KAPPA .GT. 0) FAC = -FAC
*
      FG = FAC*SQRT (1.0D 00+EPS)
      FF = FAC*SQRT (1.0D 00-EPS)
*
*   Now set up the coefficients of the confluent hypergeometric
*   functions  F (-NR+1,2*GAMMA+1;RHO)  and  F (-NR,2*GAMMA+1;RHO)
*   in the workspace arrays  TA  and  TB , respectively
*
      IF (NR .EQ. 0) THEN
         IORDR1 = 0
         IORDR2 = 0
      ELSE
         IORDR1 = NR-1
         IORDR2 = NR
      ENDIF
*
      FAC = 1.0D 00
      FACN = 1.0D 00
      A = -FNR
      AN1 = A+1.0D 00
      AN2 = A
      B = TWOGP1
      BN = B
*
      K = 0
    2 K = K+1
      FDEN = 1.0D 00/(FACN*BN)
      IF (K .LE. IORDR1) THEN
         TA(K) = AN1*FDEN
      ENDIF
      IF (K .LE. IORDR2) THEN
         TB(K) = AN2*FDEN
         A = A+1.0D 00
         AN1 = AN1*(A+1.0D 00)
         AN2 = AN2*A
         B = B+1.0D 00
         BN = BN*B
         FAC = FAC+1.0D 00
         FACN = FACN*FAC
         GOTO 2
      ENDIF
*
*   Now tabulate the function over the entire grid
*
      RG(1) = 0.0D 00
      RF(1) = 0.0D 00
      FAC = (Z+Z)/BIGN
      BIGNMK = BIGN-FKAPPA
      DO 4 I = 2,NTP
         RHO = FAC*R(I)
         RHON = RHO
         K = 0
         F1 = 1.0D 00
         F2 = 1.0D 00
    3    K = K+1
         IF (K .LE. IORDR1) THEN
            F1 = F1+TA(K)*RHON
         ENDIF
         IF (K .LE. IORDR2) THEN
            F2 = F2+TB(K)*RHON
            RHON = RHON*RHO
            GOTO 3
         ENDIF
         F1 = FNR*F1
         F2 = BIGNMK*F2
         OVLFAC = EXP (-0.5D 00*RHO)*(RHO**GAMMA)
         RG(I) = FG*OVLFAC*(F1-F2)
         RF(I) = FF*OVLFAC*(F1+F2)
    4 CONTINUE
*
*   Determine the effective maximum tabulation point based on the
*   cutoff; define the cutoff conservatively
*
      CUTOFF = ACCY*0.1D 00
*
      MTP = NTP+1
    5 MTP = MTP-1
      IF (ABS (RG(MTP)) .LT. CUTOFF) THEN
         RG(MTP) = 0.0D 00
         RF(MTP) = 0.0D 00
         GOTO 5
      ENDIF
*
      IF (MTP .EQ. NTP) WRITE (*,306) NTP,RG(NTP),CUTOFF
*
*   Compute the coefficient of R**GAMMA at the origin
*
      RG0 = FG*(FAC**GAMMA)*(FNR-BIGNMK)
*
      RETURN
*
  300 FORMAT ('DCBSRW:')
  301 FORMAT (' Principal quantum number is ',1I3)
  302 FORMAT (' Angular quantum number is 0')
  303 FORMAT (' Angular quantum number (',1I3,') is out of range for',
     :        ' principal quantum number (',1I3,')')
  304 FORMAT (' Nuclear charge (',3P,1D16.7,') is too small')
  305 FORMAT (' Nuclear charge (',3P,1D16.7,') exceeds C (',1D16.7,')')
  306 FORMAT (///' ***** Warning in SUBROUTINE DCBSRW *****'
     :         //' Radial grid of insufficient extent:'
     :          /' P(',1I4,') = ',1P,1D10.3,
     :           ', Exceeds cutoff (',1D10.3,')')
*
      END
