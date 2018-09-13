********************************************************************
*                                                                  *
      SUBROUTINE EL5(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 11  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :    N'1 = N1 (+-) 1        *
*                                           N'2 = N2 (+-) 1        *
*                                           N'3 = N3 (+-) 1        *
*                                           N'4 = N4 (+-) 1        *
*                                                                  *
*      SUBROUTINE CALLED: EL51,EL52,EL53                           *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
C      WRITE(99,101)JA,JB,JC,JD
C  101 FORMAT(5X,'    EL5    JA JB JC JD  ',4I6)
      IF(NPEEL.LE.3)RETURN
      IF(JB.LT.JC)GO TO 1
      IF(JA.GT.JD.AND.JB.GT.JD)GO TO 2
      IF(JB.GT.JC.AND.JB.LT.JD.AND.JA.LT.JC)GO TO 3
      IF(JB.GT.JC.AND.JB.GT.JD.AND.JA.GT.JC)GO TO 4
      IF(JB.GT.JC.AND.JB.GT.JD.AND.JA.LT.JC)GO TO 5
      IF(JB.GT.JC.AND.JB.LT.JD.AND.JA.GT.JC)GO TO 6
      GO TO 10
    1 CALL EL51(JJA,JJB,JA,JB,JC,JD,1,ICOLBREI)
C      WRITE(99,102)JA,JB,JC,JD
C  102 FORMAT(5X,'REZ1 go out   EL5 JA,JB,JC,JD  ',4I5)
      RETURN
    2 CALL EL51(JJA,JJB,JC,JD,JA,JB,2,ICOLBREI)
C      WRITE(99,103)JA,JB,JC,JD
C  103 FORMAT(5X,'REZ2 go out   EL5 JA,JB,JC,JD  ',4I5)
      RETURN
    3 CALL EL52(JJA,JJB,JA,JC,JB,JD,1,ICOLBREI)
C      WRITE(99,106)JA,JB,JC,JD
C  106 FORMAT(5X,'REZ1 go out   EL5 JA,JB,JC,JD  ',4I5)
      RETURN
    4 CALL EL52(JJA,JJB,JC,JA,JD,JB,2,ICOLBREI)
C      WRITE(99,105)JA,JB,JC,JD
C  105 FORMAT(5X,'REZ2 go out   EL5 JA,JB,JC,JD  ',4I5)
      RETURN
    5 CALL EL53(JJA,JJB,JA,JC,JD,JB,1,ICOLBREI)
C      WRITE(99,104)JA,JB,JC,JD
C  104 FORMAT(5X,'REZ1 go out   EL5 JA,JB,JC,JD  ',4I5)
      RETURN
    6 CALL EL53(JJA,JJB,JC,JA,JB,JD,2,ICOLBREI)
C      WRITE(99,107)JA,JB,JC,JD
C  107 FORMAT(5X,'REZ2 go out   EL5 JA,JB,JC,JD  ',4I5)
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL5 ')
      STOP
      END
