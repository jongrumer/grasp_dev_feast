************************************************************************
*                                                                      *
      SUBROUTINE SETXUV (J)
*                                                                      *
*   This  SUBROUTINE  sets  up the arrays XU and XV, for use by  the   *
*   subprograms IN and  OUT.                                           *
*                                                                      *
*   Written by Farid A Parpia, at Oxford    Last update: 17 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
CGG      PARAMETER (NNNW = 120)
*
      COMMON/DEF2/C
     :      /DEFC/DP(NNNP),DQ(NNNP)
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N
     :      /INT3/TF(NNNP),TG(NNNP),XU(NNNP),XV(NNNP)
     :      /POTE/YP(NNNP),XP(NNNP),XQ(NNNP)
     :      /SCF3/SCNSTY(NNNW),METHOD(NNNW)
*
*   Define constants
*
      DMHH = -H*0.5D 00
*
*   Set up arrays XU and XV; since XU(1), XV(1) are never used,
*   set them to some arbitrary value
*
      NM1 = N-1
      XU(1) = 0.0D 00
      XV(1) = 0.0D 00
      DO 1 I = 2,NM1
         XU(I) = DMHH*(XP(I+1)*RPOR(I+1)+XP(I)*RPOR(I))+DP(I)
         XV(I) = DMHH*(XQ(I+1)*RPOR(I+1)+XQ(I)*RPOR(I))+DQ(I)
    1 CONTINUE
!     print*,"setxuv: xu(2),xv(2)",xu(2)," ",xv(2),"myid=",myid
!     print *, "setxuv: xu(50),xv(50)",xu(50)," ", xv(50)
!     print *, "setxuv: xu(90),xv(90)",xu(90)," ", xv(90)                       
*
      XU(N) = DMHH*XP(N)*RPOR(N)+DP(N)
      XV(N) = DMHH*XP(N)*RPOR(N)+DQ(N)
*
      RETURN
      END
