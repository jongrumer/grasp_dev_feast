************************************************************************
*                                                                      *
      SUBROUTINE SPEAK (JA,JB,IA1,IB1,IA2,IB2,K,X)
*                                                                      *
*   Output MCP  coefficients and integral parameters to COMMON block   *
*   /BUFFER/. Also print these if  IBUG1 = 1 .                         *
*                                                                      *
*                                           Last Update: 14 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      CHARACTER*2 NH
*
      POINTER (PLABEL,LABEL(6,*))
      POINTER (PCOEFF,COEFF(*))
*
      COMMON/BUFFER/NBDIM,PLABEL,PCOEFF,NVCOEF
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
*
      IF (IBUG1 .NE. 0) WRITE (99,300) JA,JB,NP(IA1),NH(IA1),
     :                        NP(IB1),NH(IB1),NP(IA2),NH(IA2),
     :                        NP(IB2),NH(IB2),K,X
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
      LABEL(1,NVCOEF) = IA1
      LABEL(2,NVCOEF) = IB1
      LABEL(3,NVCOEF) = IA2
      LABEL(4,NVCOEF) = IB2
      LABEL(5,NVCOEF) = K
      COEFF(  NVCOEF) = X
*
      RETURN
*
  300 FORMAT (2(1X,1I2),4(1X,I2,A2),1X,I2,1X,1PD19.12)
*
      END
