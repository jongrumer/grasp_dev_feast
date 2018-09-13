************************************************************************
*                                                                      *
      FUNCTION SKINT (RAC,IA,IC,RBD,IB,ID,K,IW)
*                                                                      *
*   This routine evaluates transverse interaction integrals:           *
*                                                                      *
*                          (k)                                         *
*                         S   (a,c;b,d;w)                              *
*                                                                      *
*   where w = wac if IW = 1, and w= wbd if IW = 2.                     *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*               [RCI92]: ZKF.                                          *
*                                           Last update: 06 Nov 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
      LOGICAL LDBPR
      CHARACTER*2 NH
*
      DIMENSION RAC(NNNP),RBD(NNNP)
      DIMENSION TKEEP(NNNP)
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/BESS1/WIJ(2),BESSJ(2,2,NNNP),BESSN(2,2,NNNP)
     :      /DEBUGR/LDBPR(30)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /TATB/TA(NNN1),TB(NNN1),MTP
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
*
      EPSI = 1.0D-10
*
      W = WIJ(IW)
      WK = DBLE (K+K+1)
      MXRAC = MIN (MF(IA),MF(IC))
      MXRBD = MIN (MF(IB),MF(ID))
*
*           (k-1)
*  Compute Z     (rho  ; s)
*                    ac
*
      DO 1 I = 1,MXRAC
         TA(I) = RAC(I)
    1 CONTINUE
      MTP = MXRAC
      CALL ZKF (K-1,IA,IC)
*
      IF (ABS (W) .LT. EPSI) THEN
*
*   W = 0 case
*
         DO 2 I = 1,MTP
            TKEEP(I) = TB(I)
    2    CONTINUE
         CALL ZKF (K+1,IA,IC)
         MTP = MIN (MTP,MXRBD)
         TA(1) = 0.0D 00
         DO 3 I = 2,MTP
            TA(I) = RBD(I)*RPOR(I)*(TB(I)-TKEEP(I))
    3    CONTINUE
         CALL QUAD (VALU)
         SKINT = WK*VALU*0.5D 00
*
      ELSE
*
*   Finite w: see I P Grant and B J McKenzie, J Phys B: At Mol Phys,
*   13 (1980) 2671-2681
*
         DO 4 I = 1,MTP
            TKEEP(I) = TB(I)
            TA(I) = -TA(I)*BESSJ(1,IW,I)
    4    CONTINUE
         CALL ZKF (K-1,IA,IC)
*
         MTP = MIN (MTP,MXRBD)
         TA(1) = 0.0D 00
         DO 5 I = 2,MTP
            TA(I) = ((1.0D 00+BESSN(2,IW,I))
     :               *TB(I)-TKEEP(I)*BESSN(2,IW,I))*
     :                RBD(I)/(R(I)**2)*RPOR(I)
    5    CONTINUE
         CALL QUAD (VALU)
         SKINT = ((WK/W)**2) * VALU
         DO 6 I = 1,MXRBD
            TA(I) = RBD(I)*(1.0D 00+BESSJ(2,IW,I))
    6    CONTINUE
         MTP = MXRBD
         CALL ZKF (K+1,IB,ID)
         MTP = MIN (MTP,MXRAC)
         TA(1) = 0.0D 00
         DO 7 I = 2,MTP
            TA(I) = RAC(I)*(1.0D 00+BESSN(1,IW,I))
     :                    *TB(I)*R(I)*RP(I)
    7    CONTINUE
         CALL QUAD (VALU)
         SKINT = SKINT-VALU*W*W/DBLE((2*K+3)*(2*K-1))
*
      ENDIF
*
      IF (LDBPR(11)) WRITE (99,300) K,NP(IA),NH(IA),NP(IC),NH(IC),
     :                                NP(IB),NH(IB),NP(ID),NH(ID),
     :                                IW,VALU
*
      RETURN
*
  300 FORMAT ('  (',1I2,')'
     :       /'S     (',1I2,1A2,',',1I2,1A2,'|'
     :                 ,1I2,1A2,',',1I2,1A2,';',1I2,') = ',1PD19.12)
*
      END
