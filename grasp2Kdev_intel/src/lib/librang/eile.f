********************************************************************
*                                                                  *
      SUBROUTINE EILE(JA,JB,JC,JAA,JBB,JCC)
*                                                                  *
*     ------------  SECTION METWO    SUBPROGRAM 02  ------------   *
*                                                                  *
*     NO SUBROUTINE CALLED                                         *
*                                                                  *
********************************************************************
*
      JAA=JA
      JCC=JA
      IF(JAA.GT.JB)JAA=JB
      IF(JCC.LT.JB)JCC=JB
      IF(JAA.GT.JC)JAA=JC
      IF(JCC.LT.JC)JCC=JC
      IF((JA.GT.JAA).AND.(JA.LT.JCC))JBB=JA
      IF((JB.GT.JAA).AND.(JB.LT.JCC))JBB=JB
      IF((JC.GT.JAA).AND.(JC.LT.JCC))JBB=JC
      RETURN
      END
