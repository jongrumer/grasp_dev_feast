********************************************************************
*                                                                  *
      FUNCTION NMTEJJ(I2Q,I2J,J,NK,ND)
*                                                                  *
*     ------------  SECTION METWO    SUBPROGRAM 21  -------------  *
*                                                                  *
*     NO FUNCTION CALLED                                           *
*                                                                  *
********************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      DIMENSION LP(9),LG(9),LP3(27),LG3(27)
      COMMON/MTJJ/MT(63)
      COMMON/MTJJ2/MT9(6),MT11(189)
      EXTERNAL TERMJJ
      DATA LP/1,0,3,0,6,0,12,0,26/,LG/2,0,5,0,11,0,25,0,63/,
     *LP3/1,0,8,0,16,0,25,0,35,0,46,0,58,0,71,0,85,0,100,0,
     *116,0,133,0,151,0,170/,
     *LG3/7,0,15,0,24,0,34,0,45,0,57,0,70,0,84,0,99,0,115,0,
     *132,0,150,0,169,0,189/
      NMTEJJ=0
      IF(J.GT.37)RETURN
      I2V=J+1-2*I2Q
      IN=(I2Q*100+I2V)*100+I2J
      IF(J.LT.9) THEN
         JP=LP(J)
         IF(JP.EQ.0)RETURN
         JG=LG(J)
         IF(JG.EQ.0)RETURN
         JJ=JP
    2    IF(IN-MT(JJ))3,1,3
    3    JJ=JJ+1
         IF(JJ.GT.JG)RETURN
         GO TO 2
    1    NMTEJJ=JJ
      ELSEIF(J.EQ.9) THEN
        IF(MAX0(NK,ND).LT.3) THEN
          JP=1
          JG=6
          JJ=JP
    6     IF(IN-MT9(JJ))4,5,4
    4     JJ=JJ+1
          IF(JJ.GT.JG)RETURN
          GO TO 6
    5     NMTEJJ=JJ+300
        ELSE
          PRINT*, "ERROR in FUNCTION NMTEJJ"
          STOP
        END IF
      ELSE
        IL=J-10
        JP=LP3(IL)
        JG=LG3(IL)
        JJ=JP
   22   IF(IN-MT11(JJ))23,21,23
   23   JJ=JJ+1
        IF(JJ.GT.JG)RETURN
        GO TO 22
   21   NMTEJJ=JJ
      ENDIF
      RETURN
      END
