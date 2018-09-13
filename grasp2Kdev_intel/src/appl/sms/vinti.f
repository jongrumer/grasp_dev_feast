************************************************************************
*                                                                      *
      FUNCTION VINTI (J,K)
*                                                                      *
*   The value of this  function is the one-electron integral V (J,K)   *
*   for  orbitals  J, K. The analytical expression for this quantity   *
*   is given as  eq (3.23) in  F A  Parpia, M. Tong and C F Fischer,   *
*   to appear.                                                         *
*                                                                      *
*   Call(s) to: [LIB92]: DPBDT, QUAD.                                  *
*                                                                      *
*   Written by M Tong and F A Parpia,     Last revision: 15 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = 600)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      MTP = MAX (MF(J),MF(K))
*
*   Piece involving derivatives
*
      CALL DPBDT (K)
      TA(1) = 0.0D 00
      DO 1 I = 2,MTP
         TA(I) = PF(I,J)*TA(I)+QF(I,J)*TB(I)
    1 CONTINUE
      CALL QUAD (PIECE1)
      PIECE1 = PIECE1/H
*
*   Pieces not involving derivatives
*
      KPJ = NAK(J)
      KPK = NAK(K)
      IFACT1 = KPJ*(KPJ+1)-KPK*(KPK+1)
      FACT1 = 0.5D 00*DBLE (IFACT1)
      IFACT2 = -KPJ*(-KPJ+1)+KPK*(-KPK+1)
      FACT2 = 0.5D 00*DBLE (IFACT2)
      TA(1) = 0.0D 00
      DO 2 I = 2,MTP
         TA(I) = RPOR(I)*( FACT1*PF(I,J)*PF(I,K)
     :                    +FACT2*QF(I,J)*QF(I,K))
    2 CONTINUE
      CALL QUAD (PIECE2)
*
      VINTI = PIECE1-PIECE2
*
*   Debug printout
*
      IF (LDBPR(6)) WRITE (99,300) NP(J),NH(J),NP(K),NH(K),VINTI
*
      RETURN
*
  300 FORMAT (/'VINTI: V (',1I2,1A2,',',1I2,1A2,') = ',1PD19.12)
*
      END
