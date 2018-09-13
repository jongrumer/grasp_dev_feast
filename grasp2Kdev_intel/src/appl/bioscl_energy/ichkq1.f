************************************************************************
*                                                                      *
      FUNCTION ICHKQ1(JA,JB)
*                                                                      *
*   This routine is to check the occupation condition for one electron *
*   operator.                                                          *
*                                                                      *
*   Call(s) to: [LIB92]: IQ.                                           *
*                                                                      *
*   Yu Zou                                Last revision: 8/16/00       *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/DEBUGA/LDBPA(5)
     :      /ORB2/NCF,NW,PNTRIQ
*
*
      ICHKQ1=0
      K=0
      DO I=1,NW
       IQA=IQ(I,JA)
       IQB=IQ(I,JB)
       IF(IQA.NE.IQB) THEN
         K=K+1
         IF(K.GT.2) RETURN
         IF(IABS(IQA-IQB).GT.1) RETURN
       ENDIF
      ENDDO
      IF(K.EQ.2.OR.K.EQ.0) ICHKQ1=1
      RETURN
      END
