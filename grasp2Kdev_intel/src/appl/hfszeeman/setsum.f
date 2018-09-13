************************************************************************
*                                                                      *
      SUBROUTINE SETSUM(NAME,NCI,NOFFD)
*                                                                      *
*   Open the .gjhfs file on stream 111 and .ch on  29                  *
*                                                                      *
*   Call(s) to: [LIB92]: LENGTH, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
*   Updated by Per Jonsson                               28 Oct 1999   *
*   Updated by Per and Martin                                          *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      CHARACTER*256 FILNAM1,FILNAM2,FILNAM3,FILNAM4
      CHARACTER*24 NAME
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
*   File <name>.gjhfs and <name>.h or <name>.cgjhfs and 
*   <name>.ch is FORMATTED
*
      K = INDEX(NAME,' ')
      IF (NCI.EQ.0) THEN
         FILNAM1 = NAME(1:K-1)//'.ch'
         IF (NOFFD .EQ. 0) THEN
            FILNAM2 = NAME(1:K-1)//'.cgjhfs'
         ENDIF
      ELSE
         FILNAM1 = NAME(1:K-1)//'.h'
         IF (NOFFD .EQ. 0) THEN
            FILNAM2 = NAME(1:K-1)//'.gjhfs'
         ENDIF
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
      IF (NOFFD .EQ. 0) THEN
         CALL OPENFL (111,FILNAM2,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
            PRINT *, 'Error when opening',FILNAM2
            STOP
         ENDIF
      ENDIF
*
      RETURN
      END
