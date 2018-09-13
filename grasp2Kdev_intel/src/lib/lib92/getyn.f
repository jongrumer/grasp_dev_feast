************************************************************************
*                                                                      *
      FUNCTION GETYN ()
*                                                                      *
*   This  subprogram reads a response on  the default input unit; the  *
*   response must be either 'y' or 'n'. GETYN is .TRUE. if 'y' is en-  *
*   tered and .FALSE. if 'n' is entered.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 27 Aug 1992   *
*                                                                      *
************************************************************************
      COMMON/iounit/istdi,istdo,istde
      LOGICAL GETYN
      CHARACTER*1 RSPNS
*
    1 READ (istdi,'(A)') RSPNS
      IF (RSPNS .EQ. 'y') THEN
        GETYN = .TRUE.
      ELSEIF (RSPNS .EQ. 'n') THEN
        GETYN = .FALSE.
      ELSE
        WRITE(istde,*) 'Expecting <y><cr> or <n><cr> ...'
        GOTO 1
      ENDIF
*
      RETURN
      END
