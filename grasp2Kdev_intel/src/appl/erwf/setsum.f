************************************************************************
*                                                                      *
      SUBROUTINE SETSUM
*                                                                      *
*   Open the  .sum  file on stream 24.                                 *
*                                                                      *
*   Call(s) to: [LIB92]:  OPENFL.                                      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      COMMON/iounit/istdi,istdo,istde
*
*   File  erwf.sum  is FORMATTED
*
      DEFNAM = 'erwf.sum'
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
      WRITE(istde,*) 'File  erwf.sum  will be created as the '
     &             , 'ERWF SUMmary File; '
      WRITE(istde,*) 'enter another file name if this is not '
     &             , 'acceptable; null otherwise:'
      READ (*,'(A)') FILNAM
*
      IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    1 CALL OPENFL (24,FILNAM,FORM,STATUS,IERR)
      IF (IERR .NE. 0) THEN
    2    WRITE(istde,*) 'Enter a name for the ERWF SUMmary File '
     &                , 'that is to be created:'
         READ (*,'(A)') FILNAM
         IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 2
         GOTO 1
      ENDIF
*
      RETURN
      END
