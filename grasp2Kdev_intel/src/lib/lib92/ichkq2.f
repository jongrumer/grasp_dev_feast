************************************************************************
*                                                                      *
      FUNCTION ICHKQ2(JA,JB)
*                                                                      *
*   This routine is to check the occupation condition for two electron *
*   operator.                                                          *
*                                                                      *
*   Call(s) to: [LIB92]: IQ.                                           *
*                                                                      *
*   Yu Zou                                Last revision: 8/21/00       *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/DEBUGA/LDBPA(5)
     :      /ORB2/NCF,NW,PNTRIQ
*
*
      ICHKQ2=0
      K=0
      DO I=1,NW
       IQA=IQ(I,JA)
       IQB=IQ(I,JB)
       IF(IQA.NE.IQB) THEN
         K=K+1
         IF(K.GT.4) RETURN
         IF(IABS(IQA-IQB).GT.2) RETURN
       ENDIF
      ENDDO
      ICHKQ2=1
      RETURN
      END
