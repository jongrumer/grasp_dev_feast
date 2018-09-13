********************************************************************
*                                                                  *
      SUBROUTINE PERKO1(JA,BK,IK,BD,ID)
*                                                                  *
*     ------------  SECTION METWO    SUBPROGRAM 22  -------------  *
*                                                                  *
*     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
*                                               (FOR ONE SHELL)    *
*                                                                  *
*     NO SUBROUTINE CALLED                                         *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/M1/NQ1(NNNW),NQ2(NNNW)
     :      /M2/JJQ1(3,NNNW),JJQ2(3,NNNW)
      COMMON/M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
     :      /ORB4/NP(NNNW),NAK(NNNW)
      DIMENSION BK(3),IK(7),BD(3),ID(7)
      IJ=JLIST(JA)
      IK(2)=NP(IJ)
      ID(2)=IK(2)
      IK(3)=(IABS(NAK(IJ))*2)-1
      ID(3)=IK(3)
      IK(4)=NQ1(IJ)
      ID(4)=NQ2(IJ)
      IK(5)=(IK(3)+NAK(IJ)/IABS(NAK(IJ)))/2
      ID(5)=IK(5)
      IK(6)=JJQ1(3,IJ)-1
      ID(6)=JJQ2(3,IJ)-1
      IK(7)=IABS(NAK(IJ))-JJQ1(1,IJ)
      ID(7)=IABS(NAK(IJ))-JJQ2(1,IJ)
      BK(1)=HALF*DBLE(IK(7))
      BD(1)=HALF*DBLE(ID(7))
      BK(2)=HALF*DBLE(IK(6))
      BD(2)=HALF*DBLE(ID(6))
      BK(3)=-HALF*DBLE(IABS(NAK(IJ))-IK(4))
      BD(3)=-HALF*DBLE(IABS(NAK(IJ))-ID(4))
      IK(1)=NMTEJJ(IK(7),IK(6),IK(3),ID(4),IK(4))
      ID(1)=NMTEJJ(ID(7),ID(6),ID(3),ID(4),IK(4))
      RETURN
      END
