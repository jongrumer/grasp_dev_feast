************************************************************************
*                                                                      *
      FUNCTION RKINT (RAC,IA,IC,RBD,IB,ID,K,IW)
*                                                                      *
*   This routine evaluates the transverse interaction integrals.  If   *
*   IW = 0, it calulates the U(r1,r2) integral; otherwise, it calcu-   *
*   lates  R bar (k; a c | b d ; w)  with  w = wac  if  IW = 1,  and   *
*   w = wbd  if IW = 2.                                                *
*                                                                      *
*   Call(s) to: [LIB92]: QUAD.                                         *
*               [RCI92]: ZKF.                                          *
*                                           Last update: 15 Oct 1992   *
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
      MXRBD = MIN (MF(IB),MF(ID))
      MXRAC = MIN (MF(IA),MF(IC))
*
      IF (IW .EQ. 0) THEN
*
*   IW = 0
*
         DO 1 I = 1,MXRAC
            TA(I) = RAC(I)
    1    CONTINUE
         MTP = MXRAC
         CALL ZKF (K,IA,IC)
         MTP = MIN (MTP,MXRBD)
         TA(1) = 0.0D 00
         DO 2 I = 2,MTP
            TA(I) = RBD(I)*TB(I)*RPOR(I)
    2    CONTINUE
*
      ELSE
*
*   IW = 1,2
*
         DO 4 I = 1,MXRAC
            TA(I) = RAC(I)*(1.0D 00+BESSJ(1,IW,I))
    4    CONTINUE
         MTP = MXRAC
         CALL ZKF (K,IA,IC)
         MTP = MIN (MTP,MXRBD)
         TA(1) = 0.0D 00
         DO 5 I = 2,MTP
            TA(I) = RBD(I)*(1.0D 00+BESSN(1,IW,I))
     :                    *TB(I)*RPOR(I)
    5    CONTINUE
*
      ENDIF
*
      CALL QUAD (RESULT)
      RKINT = RESULT
*
*   Debug printout if option set
*
      IF (LDBPR(11)) WRITE (99,300) K,NP(IA),NH(IA),NP(IC),NH(IC),
     :                                NP(IB),NH(IB),NP(ID),NH(ID),
     :                                IW,RESULT
*
      RETURN
*
  300 FORMAT ('_ (',1I2,')'
     :       /'R     (',1I2,1A2,',',1I2,1A2,'|'
     :                 ,1I2,1A2,',',1I2,1A2,';',1I2,') = ',1PD19.12)
*
      END
