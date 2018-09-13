************************************************************************
*                                                                      *
      FUNCTION RINTDENS (I,J)
*                                                                      *
*   The value of RINTDENS is an approximation to:                      *
*                                                                      *
*                                                                      *
*      (4pi^)-1  r^-2 |P (r)*P (r) + Q (r)*Q (r) | r -> 0              *
*                     I     J       I     J                            *
*                                                                      *
*   Call(s) to: [SMS92]:  POLINT                                       *
*                                                                      *
*   Written by Per Jonsson             Last revision: 24 Dec 1992      *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION XA(3),YA(3)
      POINTER (PNTRPF,PF(NNNP,*)),(PNTRQF,QF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /DEF9/CVAC,PI
      DO 1 L = 4,2,-1
         XA(L-1) = R(L)
         YA(L-1) = (PF(L,I)*PF(L,J)+QF(L,I)*QF(L,J))/
     :             (4.0D 00*PI*R(L)*R(L))
    1 CONTINUE
      CALL POLINT(XA,YA,3,0.0D 00,DENS,DDENS)
      RINTDENS = DENS
*
      RETURN
      END
