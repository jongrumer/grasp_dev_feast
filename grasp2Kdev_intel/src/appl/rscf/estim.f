************************************************************************
*                                                                      *
      SUBROUTINE ESTIM (J)
*                                                                      *
*   This  subprogram implements  Part 1 of Algorithm 7.1 of C Froese   *
*   Fischer, Comput Phys Rep, 3 (1986) 320-321.                        *
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
*
      COMMON/DEF1/ATW,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /SCF4/EPSMIN,EPSMAX,EMIN,EMAX,ZINF
*
*   Initializations
*
      ALPHA = 1.0D 00/C
      CSQ = C*C
      NPJ = NP(J)
      NAKABS = ABS (NAK(J))
      FNREL = DBLE (NPJ-NAKABS)
      FKABS = DBLE (NAKABS)
      FKAP2 = FKABS*FKABS
*
*   Set ZINF, the asymptotic charge seen by the electron
*
*     ZINF = DBLE (IONCTY+1)
*
*   Changed on 07/06/93 by WPW
*
      ZINF = Z + DBLE (-NELEC+1)
*
*   Set the lower bound
*
      ZALPHA = ZINF*ALPHA
      IF (ZALPHA .LT. FKABS) THEN
         GAMMA = SQRT (FKAP2-ZALPHA*ZALPHA)
         EBYM = 1.0D 00
     :          /SQRT ( 1.0D 00
     :               +(ZALPHA/(GAMMA+FNREL+0.5D 00))**2)
         EPSMIN = (1.0D 00-EBYM)*CSQ
      ELSE
         EPSMIN = 0.25D 00*CSQ/DBLE (NPJ*NPJ)
      ENDIF
      EMIN = EPSMIN
*
*   Set the upper bound
*
      ZALPHA = Z*ALPHA
      IF (ZALPHA .LT. FKABS) THEN
         GAMMA = SQRT (FKAP2-ZALPHA*ZALPHA)
         EBYM = 1.0D 00
     :          /SQRT ( 1.0D 00
     :                 +(ZALPHA/(GAMMA+FNREL-0.5D 00))**2)
         EPSMAX = (1.0D 00-EBYM)*CSQ
      ELSE
         EPSMAX = CSQ+CSQ
      ENDIF
      EMAX = EPSMAX
*
      RETURN
      END
