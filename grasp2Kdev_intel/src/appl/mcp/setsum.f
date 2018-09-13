************************************************************************
*                                                                      *
      SUBROUTINE SETSUM (fullname)
*                                                                      *
*   Open the  .sum  file on stream 24.                                 *
*                                                                      *
*   Call(s) to: [LIB92]:  OPENFL.                                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 11 Nov 1992   *
*   Modified by Xinghong He               Last revision:  3 Jul 1998   *
*
*  File shared by mcpblk, mcpmpi
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H,O-Z)
      CHARACTER*(*) fullname
      CHARACTER(LEN = LEN (fullname)) FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS

      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      WRITE (istde,*) 'File ', fullname,' will be created as the'
     &,        ' GENMCP SUMmary File;'
      WRITE (istde,*) 'enter another file name if this is not '
     &,        'acceptable; null otherwise:'
      READ (*,'(A)') FILNAM
*
      IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = fullname
*
    1 CALL OPENFL (24, FILNAM, FORM, STATUS, IERR)
      IF (IERR .NE. 0) THEN
    2    WRITE (istde,*) 'Enter a name for the GENMCP SUMmary'
     &,           ' File that is to be created:'
         READ (*,'(A)') FILNAM
         IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 2
         GOTO 1
      ENDIF
*
      RETURN
      END
