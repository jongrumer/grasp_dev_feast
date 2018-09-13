************************************************************************
*                                                                      *
      SUBROUTINE TAIL (IORB,P,Q,JP,MTP)
*                                                                      *
*   This subroutine begins the inward integration of the homogeneous   *
*   Dirac radial equation. With only minor modifications, the series   *
*   given by  J E Sienkiewicz  and W E Baylis, J Phys B: At Mol Phys   *
*   20 (1987) 5145-5156, p 5155, is used.                              *
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
*
      DIMENSION P(NNNP),Q(NNNP)
*
      COMMON/DEF0/TENMAX,EXPMAX,EPSI,PRECIS
     :      /DEF2/C
     :      /DEF4/ACCY,NSCF,NSIC,NSOLV
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
*
*   Initialize
*
      EPS = ACCY*0.1D 00
      BIGE = -E(IORB)
      BIGEBC = BIGE/C
      BGEBC2 = BIGEBC/C
      BETA = SQRT (-BIGE*(2.0D 00+BGEBC2))
      FK = DBLE (NAK(IORB))
*
*   Find MTP
*
      I = JP
      NM4 = N-4
      RJP = R(JP)
      QEM = 0.25D 00*EXPMAX
    1 I = I+1
      IF (I .LE. NM4) THEN
         IF (BETA*(R(I)-RJP) .GT. QEM) THEN
            MTP = I
         ELSE
            GOTO 1
         ENDIF
      ELSE
         WRITE (*,300)
         MTP = NM4
      ENDIF
*
*   Compute offset for exponential function
*
      OFFSET = BETA*R(MTP)
*
*   Tabulate tail points
*
      DO 3 I = MTP,N
*
         T = 2.0D 00*BETA*R(I)
*
         ZLOC = YP(I)
         GAM2 = FK**2-(ZLOC/C)**2
         ZLBB = ZLOC/BETA
         DNU = ZLBB*(1.0D 00+BGEBC2)
         FKMZBB = FK-ZLBB
*
         M = -1
         SUMP = 0.0D 00
         SUMQ = 0.0D 00
    2    M = M+1
         EM = DBLE (M)
         IF (M .EQ. 0) THEN
            CM = 1.0D 00
            EMFACT = 1.0D 00
         ELSE
            CM = CM*(GAM2-(DNU-EM)**2)
            EMFACT = EMFACT*EM
         ENDIF
         OVLTRM = CM*((T**(DNU-EM))/EMFACT)
         PTERM = OVLTRM*(FKMZBB+EM)*BETA
         QTERM = OVLTRM*(FKMZBB-EM)*BIGEBC
         SUMP = SUMP+PTERM
         SUMQ = SUMQ+QTERM
         IF ( (ABS (PTERM/SUMP) .GE. EPS) .OR.
     :        (ABS (QTERM/SUMQ) .GE. EPS) ) GOTO 2
         EXPTRM = EXP (-0.5D 00*T+OFFSET)
         P(I) = SUMP*EXPTRM
         Q(I) = SUMQ*EXPTRM
         IF (P(I) .EQ. 0.0D 00) THEN
            LOC = I+1
            GOTO 4
         ENDIF
    3 CONTINUE
      LOC = N+1
*
    4 DO 5 I = LOC,N
         P(I) = 0.0D 00
         Q(I) = 0.0D 00
    5 CONTINUE
*
      RETURN
*
  300 FORMAT ('TAIL: Grid may be of insufficient extent')
*
      END
