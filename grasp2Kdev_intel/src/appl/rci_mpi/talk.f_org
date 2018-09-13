************************************************************************
*                                                                      *
      SUBROUTINE TALK (JA,JB,NU,IA,IB,IC,ID,ITYPE,COEF)
*                                                                      *
*   Print  coefficients  and  integral  parameters  if IBUG1 > 0 and   *
*   write to disk.                                                     *
*                                                                      *
*                                           Last update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      CHARACTER*2 NH
*
      POINTER (PLABEL,LABEL(6,1))
      POINTER (PCOEFF,COEFF(1))
*
      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
*
CGG      PARAMETER (KEYORB = 121)
      PARAMETER (KEY = KEYORB)
*
*   Print coefficient if requested
*
      IF (IBUG1 .NE. 0) WRITE (99,300) JA,JB,
     :                  NP(IA),NH(IA),NP(IB),NH(IB),
     :                  NP(IC),NH(IC),NP(ID),NH(ID),
     :                  NU,ITYPE,COEF
*
*   Increment counter
*
      NVCOEF = NVCOEF+1
*
*   Ensure that arrays are of adequate size; reallocate if necessary
*
      IF (NVCOEF .GT. NBDIM) CALL ALCBUF (2)
*
*   Store integral indices and coefficient in COMMON/BUFFER/
*
      LABEL(1,NVCOEF) = IA
      LABEL(2,NVCOEF) = IB
      LABEL(3,NVCOEF) = IC
      LABEL(4,NVCOEF) = ID
      LABEL(5,NVCOEF) = NU
      LABEL(6,NVCOEF) = ITYPE
      COEFF(  NVCOEF) = COEF
*
      RETURN
*
  300 FORMAT (2(1X,1I2),4(1X,I2,A2),1X,1I2,1X,1I2,1X,1PD19.12)
*
      END
