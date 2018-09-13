************************************************************************
*                                                                      *
      FUNCTION ESTRMS (APARM,CPARM)
*                                                                      *
*   Determines the root mean square radius for a Fermi nucleus given   *
*   the parameters `c' (CPARM) and `a' (APARM). We use the formalism   *
*   developed in F. A. Parpia and A. K. Mohanty ``Relativistic basis   *
*   set calculations for atoms with Fermi nuclei'' Phys Rev A (1992)   *
*   in press.                                                          *
*                                                                      *
*   Call(s) to: SKFUN.                                                 *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
*
      PI = 4.0D 00*ATAN (1.0D 00)
      SQTBF = SQRT (3.0D 00/5.0D 00)
*
      ABC = APARM/CPARM
      PABC = PI*ABC
      CBAM = -CPARM/APARM
      DNUMER =  1.0D 00
     :         +(10.0D 00/3.0D 00)
     :          *PABC**2
     :         +(7.0D 00/3.0D 00)
     :          *PABC**4
     :         -120.0D 00
     :          *ABC**5
     :           *SKFUN (5,CBAM)
      DDENOM =  1.0D 00
     :          +PABC**2
     :          -6.0D 00
     :           *ABC**3
     :            *SKFUN (3,CBAM)
      ESTRMS = CPARM*SQTBF*SQRT (DNUMER/DDENOM)
*
      RETURN
      END
