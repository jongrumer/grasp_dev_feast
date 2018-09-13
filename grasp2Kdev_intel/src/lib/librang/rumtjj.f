********************************************************************
*                                                                  *
      SUBROUTINE RUMTJJ(KNT,JJ,LQ,LV,L)
*                                                                  *
*   ---------------  SECTION SQJJ  SUBPROGRAM 16  --------------   *
*                                                                  *
*     NO SUBROUTINE CALLED                                         *
*                                                                  *
********************************************************************
*
      COMMON/MTJJ/MT(63)
      COMMON/MTJJ2/MT9(6),MT11(189)
      IF(JJ.LT.9) THEN
        KT=MT(KNT)
      ELSE IF(JJ.EQ.9) THEN
        IF(KNT.GT.300) THEN
          KNTMIN=KNT-300
          KT=MT9(KNTMIN)
        ELSE
          PRINT*, "ERROR in RUMTJJ"
          STOP
CGG          CALL RUMT67(KNT,NR,LQ,LS,L)
CGG          RETURN
        ENDIF
      ELSE
        KT=MT11(KNT)
      ENDIF
      LQ=JTHN(KT,3,100)
      LV=JTHN(KT,2,100)
      L=JTHN(KT,1,100)
      RETURN
      END
