********************************************************************
*                                                                  *
      SUBROUTINE EL4(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 09  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :      N'1 = N1 +- 1        *
*                                             N'2 = N2 +- 1        *
*                                             N'3 = N3 -+ 2        *
*                                                                  *
*     SUBROUTINE CALLED: EL41                                      *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
      IF(NPEEL.LE.2)RETURN
      IF(JA.EQ.JB)GO TO 1
      IF(JC.EQ.JD)GO TO 2
      GO TO 10
    1 CALL EL41(JJA,JJB,JC,JD,JA,1,JA,JB,JC,JD,ICOLBREI)
      RETURN
    2 CALL EL41(JJA,JJB,JA,JB,JC,2,JA,JB,JC,JD,ICOLBREI)
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL4 ')
      STOP
      END
