********************************************************************
*                                                                  *
      FUNCTION ITJJ2(IK,ID,KG1,BK,BD,IBT,BT,ITP,ITG,IQ)
*                                                                  *
*   ---------------  SECTION SQJJ  SUBPROGRAM 07  --------------   *
*                                                                  *
*     FUNCTION CALLED: ITTK                                        *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      COMMON/GLCON/ZERO,TENTH,HALF,ONE,TWO,THREE,FOUR,SEVEN,ELEVEN,EPS
      COMMON/RIBOJJ/IMPTJJ(63),IMGTJJ(63),IMPNJJ(63),IMGNJJ(63)
      COMMON/RIBOJJ9/IMPTJJ9(6),IMGTJJ9(6),IMPNJJ9(6),IMGNJJ9(6)
      COMMON/RIBOJJ11/
     :   IMPTJJ11(189),IMGTJJ11(189),IMPNJJ11(189),IMGNJJ11(189)
      DIMENSION ID(7),IK(7),IBT(7),BT(3),BD(3),BK(3)
      ITJJ2=0
      IF(ID(3).GT.37) RETURN
      IF(ITTK(ID(6),IK(6),KG1).EQ.0)RETURN
      ITK=IK(1)
      ITD=ID(1)
      IF(ID(3).LT.9) THEN
        ITP1=IMPNJJ(ITK)
        ITP=IMPTJJ(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNJJ(ITK)
        ITG=IMGTJJ(ITD)
      ELSEIF(ID(3).EQ.9) THEN
        IF(ITK.GT.300) THEN
          IF(ITD.LT.300) CALL MES(52)
          IF(ID(4).GT.2) CALL MES(12)
          IF(IK(4).GT.2) CALL MES(12)
          ITK=ITK-300
          ITD=ITD-300
          ITP1=IMPNJJ9(ITK)
          ITP=IMPTJJ9(ITD)
          IF(ITP1.NE.ITP)RETURN
          ITG1=IMGNJJ9(ITK)
          ITG=IMGTJJ9(ITD)
        ELSE
          PRINT*, "ERROR in ITJJ2"
          STOP
        END IF
      ELSE
        IF(ID(4).GT.2) CALL MES(12)
        IF(IK(4).GT.2) CALL MES(12)
        ITP1=IMPNJJ11(ITK)
        ITP=IMPTJJ11(ITD)
        IF(ITP1.NE.ITP)RETURN
        ITG1=IMGNJJ11(ITK)
        ITG=IMGTJJ11(ITD)
      ENDIF
      IF(ITG1.NE.ITG)RETURN
      ITJJ2=1
      IBT(2)=ID(2)
      IBT(3)=ID(3)
      IBT(4)=ID(4)+IQ
      BT(3)=BD(3)+HALF*DBLE(IQ)
      RETURN
      END
