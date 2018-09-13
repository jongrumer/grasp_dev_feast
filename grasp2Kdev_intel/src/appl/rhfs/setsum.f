************************************************************************
*                                                                      *
      SUBROUTINE SETSUM(NAME,NCI)
*                                                                      *
*   Open the  .sum  files on stream 24 and 29                          *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
*   Updated by Per Jonsson                               28 Oct 1999   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM1,FILNAM2
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
*   File  hfs92.sum  is FORMATTED
*
      K = INDEX(NAME,' ')
      IF (NCI.EQ.0) THEN
         FILNAM1 = NAME(1:K-1)//'.ch'
         FILNAM2 = NAME(1:K-1)//'.choffd'
      ELSE
         FILNAM1 = NAME(1:K-1)//'.h'
         FILNAM2 = NAME(1:K-1)//'.hoffd'
      ENDIF
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      CALL OPENFL (29,FILNAM1,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening',FILNAM1
         STOP
      ENDIF
*
      CALL OPENFL (24,FILNAM2,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Error when opening',FILNAM2
         STOP
      ENDIF
*
      RETURN
      END
