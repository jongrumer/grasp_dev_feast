************************************************************************
*                                                                      *
      SUBROUTINE SETSUM(NAME,NCI)
*                                                                      *
*   Open the  .sum  file on stream 24.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
*   File  ris.sum  is FORMATTED
*
      K = INDEX(NAME,' ')
      IF (NCI.EQ.0) THEN
         FILNAM = NAME(1:K-1)//'.ci'
      ELSE
         FILNAM = NAME(1:K-1)//'.i'
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
    1 CALL OPENFL (24,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening',FILNAM
         STOP
      ENDIF
*
      RETURN
      END
