********************************************************************
*                                                                  *
      SUBROUTINE REC3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,REC)
*                                                                  *
*   ---------------  SECTION REC    SUBPROGRAM 07  --------------  *
*                                                                  *
*     SUBROUTINE CALLED:  RECO3                                    *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      IF((JA3.GT.JA1).AND.(JA3.GT.JA2))GO TO 1
      IF((JA3.LT.JA1).AND.(JA3.LT.JA2))GO TO 2
      IF(JA1-JA2)23,10,33
   23 CALL RECO3(JA1,JA3,JA2,K1,K3,K2,IRE,IAT,REC)
      IFAZ=K1-K2-K3
      IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
      RETURN
   33 CALL RECO3(JA2,JA3,JA1,K2,K3,K1,IRE,IAT,REC)
      IF((K1/2)*2.NE.K1)REC=-REC
      RETURN
    1 IF(JA1-JA2)21,10,31
   21 CALL RECO3(JA1,JA2,JA3,K1,K2,K3,IRE,IAT,REC)
      RETURN
   31 CALL RECO3(JA2,JA1,JA3,K2,K1,K3,IRE,IAT,REC)
      IFAZ=K1+K2-K3
      IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
      RETURN
    2 IF(JA1-JA2)22,10,32
   22 CALL RECO3(JA3,JA1,JA2,K3,K1,K2,IRE,IAT,REC)
      IF((K3/2)*2.NE.K3)REC=-REC
      RETURN
   32 CALL RECO3(JA3,JA2,JA1,K3,K2,K1,IRE,IAT,REC)
      IFAZ=K1+K2+K3
      IF((IFAZ/4)*4.NE.IFAZ)REC=-REC
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN REC')
      STOP
      END
