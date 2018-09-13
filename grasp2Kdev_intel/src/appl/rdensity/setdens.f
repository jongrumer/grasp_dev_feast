************************************************************************
*                                                                      *
      SUBROUTINE SETDENS(NAME,NCI)
*                                                                      *
*   Open the  .d    file on stream 35.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by J. Ekman                                  23 Nov 2013   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
*   File  sms92.sum  is FORMATTED
*
      K = INDEX(NAME,' ')
      IF (NCI.EQ.0) THEN
         FILNAM = NAME(1:K-1)//'.cd'
      ELSE
         FILNAM = NAME(1:K-1)//'.d'
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
    1 CALL OPENFL (35,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening',FILNAM
         STOP
      ENDIF
*
      RETURN
      END
