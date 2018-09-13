************************************************************************
*                                                                      *
      FUNCTION LDIGIT (CST)
*                                                                      *
*   .TRUE.  if  CST  is the ASCII representation of a decimal digit;   *
*   .FALSE. otherwise.                                                 *
*                                                                      *
*   Written by Farid A. Parpia             Last revised: 16 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL LDIGIT
      CHARACTER*1 CDGT,CST
*
      DIMENSION CDGT(0:9)
*
      DATA CDGT /'0','1','2','3','4','5','6','7','8','9'/
*
      DO 1 I = 0,9
         IF (CST .EQ. CDGT(I)) THEN
            LDIGIT = .TRUE.
            GOTO 2
         ENDIF
    1 CONTINUE
      LDIGIT = .FALSE.
*
    2 RETURN
      END
