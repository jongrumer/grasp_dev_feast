************************************************************************
*                                                                      *
      SUBROUTINE SUMMRY (NUNIT)
*                                                                      *
*   Prints  a summary of the complete  list of subshell radial wave-   *
*   functions on NUNIT.                                                *
*                                                                      *
*   Call(s) to: [LIB92]:                                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 SOURCE
      CHARACTER*2 NH
*
      POINTER (PNTRPF,PF(NNNP,*))
      POINTER (PNTRQF,QF(NNNP,*))
*
      COMMON/GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /ORB1/E(NNNW),GAMA(NNNW)
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
     :      /WAVE/PZ(NNNW),PNTRPF,PNTRQF,MF(NNNW)
     :      /WHFROM/SOURCE(NNNW)
*
      WRITE (NUNIT,300)
*
      
      DO 1 I = 1,NW
         LENTH = LEN_TRIM (SOURCE(I))
         WRITE (NUNIT,301) NP(I),NH(I),E(I),PZ(I),GAMA(I),PF(2,I)
     :         ,QF(2,I),MF(I), SOURCE(I)(1:LENTH)
C         WRITE (NUNIT,302) SOURCE(I)(1:LENTH)
    1 CONTINUE
*
      RETURN
*
  300 FORMAT ('Shell',6x,'e',11X,'p0',8X,
     :        'gamma',8X,'P(2)',7X,'Q(2)',6X,'MTP', '  SRC'/)
  301 FORMAT (1X,I2,A2,5D12.4,I5,2x,A3)
C  302 FORMAT ('      Source: ',A)
*
      END
