************************************************************************
      SUBROUTINE SETSUM (filnam)
      IMPLICIT NONE
************************************************************************
      CHARACTER(LEN=*) filnam
      CHARACTER(LEN=*), PARAMETER:: form = 'FORMATTED', STATUS = 'NEW'
      INTEGER istdi, istdo, istde,  ierr
      COMMON/iounit/istdi,istdo,istde
        
      CALL openfl (24, filnam, form, status, ierr)
      IF (ierr .NE. 0) THEN
         WRITE (istde,*) 'Error when opening ', filnam
         STOP
      ENDIF

      RETURN
      END
