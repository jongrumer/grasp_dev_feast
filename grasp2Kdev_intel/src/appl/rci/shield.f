************************************************************************
*                                                                      *
      FUNCTION SHIELD (J)
*                                                                      *
*   This  routine estimates the screening (or shielding) for orbital   *
*   J according to the relativistic hydrogenic energy formula.         *
*                                                                      *
*   Written by Yu Zou, at Vanderbilt University: 02 Mar 2000           *
*                                                                      *
*   Formula:                                                           *
*    epsilon = E(n,k)/Z^2_{eff}                                        *
*    beta = alpha Z_{eff}                                              *
*    epsilon = 1/beta^{2} [(1+beta^2/nu^2)^{-1/2} - 1]                 *
*    nu = n + (k^2 - beta^2)^{1/2} - |k|                               *
*   where E(n,k) is the orbital energy for the specific principal      *
*   quantum n and k quantum number, alpha is the fine structure        *
*   constant.                                                          *
*   Note that SHIELD is set to -C if Zeff is out of range. This will   *
*   ensure that the self-energy for such orbital is zero.              *
*   Reference: M J Seaton, Rep. Prog. Phys., Vol.46, P167,1983         *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF2/C
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
c
c a=[(1+\alpha^2 E)^{-2} - 1]^{-1/2}
c note E(J) = -E(n,k)
      a=1-E(J)/(C*C)
      a=1/(a*a)-1
      if(a.lt.0.0) then
        shield=-c
        return
      endif
      a=1/sqrt(a)
c
      aa=1+a*a
      bb=a*dble(np(j)-iabs(nak(j)))
      cc=dble((np(j)-iabs(nak(j)))**2-nak(j)**2)
      t=bb*bb-aa*cc
      if(t.lt.0.0) then
        shield=-c
        return
      endif
      beta=bb+sqrt(t)
      beta=beta/aa
      shield =z-beta*c
c     print *, j, '  Initial Screen Factor = ',shield
      return
      END
