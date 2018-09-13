********************************************************************
*                                                                  *
      SUBROUTINE EL3(JJA,JJB,JA,JB,JC,JD,ICOLBREI)
*                                                                  *
*   --------------  SECTION METWO    SUBPROGRAM 05  -------------  *
*                                                                  *
*     THIS PACKAGE DETERMINES THE VALUES OF MATRIX ELEMENTS        *
*     OF TWO PARTICLE OPERATOR IN CASE :       N'1 = N1 - 1        *
*                                              N'2 = N2 + 1        *
*                                                                  *
*     SUBROUTINE CALLED: EL31,EL32,EL33                            *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
      COMMON /M3/JLIST(NNNW),KLIST(NNNW),NPEEL,NCORE
      COMMON/TRK/BD1(3),BD2(3),BK1(3),BK2(3),
     *ID1(7),ID2(7),IK1(7),IK2(7)
C      WRITE(99,101)
C  101 FORMAT(5X,'    EL3')
      IF(NPEEL.LE.1)RETURN
      IF(JB.EQ.JD)GO TO 1
      IF(JA.EQ.JC)GO TO 2
      IF(JA.EQ.JD)GO TO 3
      IF(JB.EQ.JC)GO TO 4
      GO TO 10
    1 IF(JA.EQ.JB.OR.JC.EQ.JB)GO TO 11
      CALL EL33(JJA,JJB,JC,JA,JB,1,JA,JB,JC,JD,ICOLBREI)
      RETURN
    2 IF(JB.EQ.JA.OR.JD.EQ.JA)GO TO 12
      CALL EL33(JJA,JJB,JD,JB,JA,1,JA,JB,JC,JD,ICOLBREI)
      RETURN
    3 IF(JB.EQ.JA.OR.JC.EQ.JA)GO TO 13
      CALL EL33(JJA,JJB,JC,JB,JA,2,JA,JB,JD,JC,ICOLBREI)
      RETURN
    4 IF(JA.EQ.JB.OR.JD.EQ.JB)GO TO 14
      CALL EL33(JJA,JJB,JD,JA,JB,2,JA,JB,JD,JC,ICOLBREI)
      RETURN
C
C     JA JC differents  and JB JD the same
C
   11 IF(JA.EQ.JC)GO TO 10
      IF(JC.NE.JB)GO TO 21
      CALL EL31(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
      RETURN
   21 CALL EL32(JJA,JJB,JC,JA,JA,JB,JC,JD,ICOLBREI)
      RETURN
C
C     JA JC differents  and JB JD the same
C
   12 IF(JB.EQ.JD)GO TO 10
      IF(JD.NE.JA)GO TO 22
      CALL EL31(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
      RETURN
   22 CALL EL32(JJA,JJB,JD,JB,JA,JB,JC,JD,ICOLBREI)
      RETURN
C
C     JA JC differents  and JB JD the same
C
   13 IF(JB.EQ.JC)GO TO 10
      IF(JC.NE.JD)GO TO 23
      CALL EL31(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
      RETURN
   23 CALL EL32(JJA,JJB,JC,JB,JA,JB,JC,JD,ICOLBREI)
      RETURN
C
C     JA JC differents  and JB JD the same
C
   14 IF(JA.EQ.JD)GO TO 10
      IF(JD.NE.JB)GO TO 24
      CALL EL31(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
      RETURN
   24 CALL EL32(JJA,JJB,JD,JA,JA,JB,JC,JD,ICOLBREI)
      RETURN
   10 WRITE(99,100)
  100 FORMAT(5X,'ERRO IN EL3  PMGG RAGG')
      STOP
      END
