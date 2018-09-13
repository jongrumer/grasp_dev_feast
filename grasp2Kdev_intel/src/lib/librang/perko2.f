********************************************************************
*                                                                  *
      SUBROUTINE PERKO2(JA1,JA2,JA3,JA4,I)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 23  -------------  *
*                                                                  *
*     INTERFACE BETWEEN "GRASP" AND BOLCK "SQ"                     *
*                                               (GENERAL CASE)     *
*                                                                  *
*     SUBROUTINE CALLED: PERKO1                                    *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
      COMMON/TRK2/BD3(3),BD4(3),BK3(3),BK4(3),
     *ID3(7),ID4(7),IK3(7),IK4(7)
      EXTERNAL GLCONS,RIBJJ
      CALL PERKO1(JA1,BK1,IK1,BD1,ID1)
      IF(I.EQ.1)RETURN
      CALL PERKO1(JA2,BK2,IK2,BD2,ID2)
      IF(I.EQ.2)RETURN
      CALL PERKO1(JA3,BK3,IK3,BD3,ID3)
      IF(I.EQ.3)RETURN
      CALL PERKO1(JA4,BK4,IK4,BD4,ID4)
      RETURN
      END
