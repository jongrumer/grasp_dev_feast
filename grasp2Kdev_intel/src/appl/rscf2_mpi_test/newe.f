************************************************************************
*                                                                      *
      SUBROUTINE NEWE (J,SGN,NPRIME,MX,DELEPS,FAIL,INV)
*                                                                      *
*   This  subroutine implements Part 2 of Algorithm 7.1 in  C Froese   *
*   Fischer,  Comput  Phys  Rep, 3 (1986) 273-326. (The present code   *
*   actually  implements  the version used in her program MCHF where   *
*   differences occur.)                                                *
*                                                                      *
*   Call(s) to: [RSCF92]: OUTBND.                                      *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last Update: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL FAIL,OUTBND
*
      COMMON/INT2/P0,Q0,P(NNNP),Q(NNNP),MTP0
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /SCF4/EPSMIN,EPSMAX,EMIN,EMAX,ZINF
*
      PARAMETER (P02    = 2.0D-02,
     :           P05    = 5.0D-02,
     :           P001   = 1.0D-03,
     :           P00001 = 1.0D-05,
     :           D2P5   = 2.5D 00)
*
*   Determine if the iterative process has succeeded
*
      IF ((SGN .GT. 0.0D 00) .AND. (MX .EQ. 0)) THEN
         FAIL = .FALSE.
      ELSE
         FAIL = .TRUE.
      ENDIF
*
*   Inversion counter
*
      INV = 0
*
*   If unsuccessful, obtain a new estimate of the eigenvalue
*
      IF (FAIL) THEN
*
*   Define quantities used later
*
         EPS = E(J)
         ABDELE = ABS (DELEPS)
         DELEBE = ABS (DELEPS/EPS)
         DN = DBLE (NP(J))
         DNPRME = DBLE (NPRIME)
         DNPN25 = (DNPRME/DN)**D2P5
*
         IF ( ( (ABS (MX) .EQ. 1) .AND. (DELEBE .GT. P02   ) ) .OR.
     :        ( (     MX  .EQ. 0) .AND. (DELEBE .GE. P00001)
     :                            .AND. (ABDELE .GE. P001  ) ) ) THEN
    1                                            ETRY = EPS+DELEPS
            IF (OUTBND (ETRY) .AND. (MX .NE. 0)) ETRY = EPS*DNPN25
            IF (OUTBND (ETRY))                   ETRY = EPS-DELEPS
            IF (OUTBND (ETRY)) THEN
               DELEPS = 0.5D 00*DELEPS
               GOTO 1
            ENDIF
         ELSEIF (MX .EQ. 0) THEN
            ETRY = EPS
            INV = 1
            P0 = -P0
            DO 2 I = 1,MTP0
               P(I) = -P(I)
               Q(I) = -Q(I)
    2       CONTINUE
            FAIL = .FALSE.
         ELSEIF (SGN .LT. 0.0D 00) THEN
            ETRY = EPS*DNPN25
            IF (OUTBND (ETRY)) ETRY = 0.5D 00*(EPSMAX+EPS)
            IF (OUTBND (ETRY)) ETRY = 0.5D 00*(EPSMIN+EPS)
         ELSEIF (MX .LT. 0) THEN
            DELTA = 1.0D 00-EPS/EPSMAX
            EPSMAX = EPS
            IF (DELTA .LT. P05) THEN
               EMAX = EMAX*DNPN25
            ELSE
               EMAX = EPS*DNPN25
            ENDIF
            ETRY = EMAX
            IF (OUTBND (ETRY)) ETRY = 0.5D 00*(EPSMAX+EPSMIN)
         ELSEIF (MX .GT. 0) THEN
            DELTA = 1.0D 00-EPSMIN/EPS
            EPSMIN = EPS
            IF (DELTA .LT. P05) THEN
               EMIN = EMIN*DNPN25
            ELSE
               EMIN = EPS*DNPN25
            ENDIF
            ETRY = EMIN
            IF (OUTBND (ETRY)) ETRY = 0.5D 00*(EPSMAX+EPSMIN)
         ENDIF
         E(J) = ETRY
      ENDIF
*
      RETURN
      END
