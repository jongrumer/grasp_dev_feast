************************************************************************
*                                                                      *
      SUBROUTINE OUT (J,JP,P,Q)
*                                                                      *
*   This subroutine carries out the step-by-step outward integration   *
*   of a pair of inhomogeneous Dirac radial equations.                 *
*                                                                      *
*   arguments:                                                         *
*                                                                      *
*      J:   (Input) Orbital index of function to be computed           *
*      JP:  (Input) The join point; the outward integration stops      *
*           at this tabulation index                                   *
*      P,Q: (Input and output) on input, elements 1 to 3 of both       *
*           arrays must be tabulated; on output, the arrays are        *
*           tabulated up to point JP                                   *
*                                                                      *
*   Written by Farid A Parpia, at Oxford   Last updated: 08 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION P(NNNP),Q(NNNP)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
     :      /ORB4/NP(NNNW),NAK(NNNW)
*
*   One global initialization
*
      DKHB2 = 0.5D 00*H*DBLE (NAK(J))
*
*   Tabulate P(r) and Q(r) by step-by-step integration
*
*   Initializations: set quantities for I = 3
*
      I = 3
      DKHB2F = DKHB2*RPOR(I)
      CPI = 1.0D 00+DKHB2F
      CMI = 1.0D 00-DKHB2F
      PI = P(I)
      QI = Q(I)
      TFI = TF(I)
      TGI = TG(I)
*
*   March out to from I = 4 to I = JP
*
!XHH Use doo-loop
!    1 I = I+1
!      IF (I .LE. JP) THEN
      DO i = 4, jp
         CPIM1 = CPI
         CMIM1 = CMI
         DKHB2F = DKHB2*RPOR(I)
         CPI = 1.0D 00+DKHB2F
         CMI = 1.0D 00-DKHB2F
         PIM1 = PI
         QIM1 = QI
         TFIM1 = TFI
         TGIM1 = TGI
         TFI = TF(I)
         TGI = TG(I)
         UCIM1 = CMIM1*PIM1-TFIM1*QIM1+XU(I-1)
         UDIM1 = CPIM1*QIM1-TGIM1*PIM1+XV(I-1)
         UEI = CPI*CMI-TFI*TGI
         PI = (CMI*UCIM1-TFI*UDIM1)/UEI
         QI = (CPI*UDIM1-TGI*UCIM1)/UEI
         P(I) = PI
         Q(I) = QI
!         GOTO 1
!      ENDIF
      ENDDO
*
      RETURN
      END

