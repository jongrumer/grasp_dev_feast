************************************************************************
*                                                                      *
      FUNCTION LENGTH (STRING)
*                                                                      *
*   Determines the location of the  rightmost character that is nei-   *
*   ther a blank nor a null.                                           *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Sep 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*(*) STRING
*
      LRIGHT = LEN (STRING)
*
      DO 1 I = LRIGHT,1,-1
         IF (STRING(I:I) .NE. ' ') THEN
            LENGTH = I
            GOTO 2
         ENDIF
    1 CONTINUE
      LENGTH = 0
*
    2 RETURN
      END
