************************************************************************
*                                                                      *
      SUBROUTINE VACPOL
*                                                                      *
*   This  routine controls the setting up of the vacuum polarization   *
*   potential for the given nuclear charge distribution at each grid   *
*   point using the analytic functions defined by  L Wayne Fullerton   *
*   and G A Rinker Jr in Phys Rev A 13 (1976) 1283. The potential is   *
*   accumulated  in  array TB(I), I = 1, ..., N  which  is in COMMON   *
*   block /TATB/ .                                                     *
*                                                                      *
*   Call(s) to: [RCI92]: VAC2, VAC4.                                   *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 09 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /NPAR/PARM(2),NPARM
     :      /NCDIST/ZDIST(NNNP)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
*
*   Redefine  ZDIST  to be  rho*r*r'
*
      DO 1 I = 1,MTP
         ZDIST(I) = ZDIST(I)*R(I)*RP(I)
    1 CONTINUE
*
*   Second-order vacuum polarisation potential; returned in
*   array TB
*
      CALL VAC2
*
*   Fourth-order vacuum polarization potential; returned in
*   array TA
*
      CALL VAC4
*
*   If option 7 is set, use user-defined vacuum polarization
*   potential
*
      RETURN
      END
