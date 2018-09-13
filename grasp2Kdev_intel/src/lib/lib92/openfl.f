************************************************************************
*                                                                      *
      SUBROUTINE OPENFL (NFILE,FILNAM,RFORM,RSTAT,IERR)
*                                                                      *
*   Issues OPEN for file with unit number NFILE, name FILNAM, format   *
*   RFORM, status RSTAT.  If this is successful the head is position-  *
*   ed to the beginning of the file and IERR is 0; otherwise IERR is   *
*   set to 1.                                                          *
*                                                                      *
*   Call(s) to: [LIB92]:  none.                                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 05 Oct 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*(*) FILNAM, RFORM, RSTAT
      COMMON/iounit/istdi,istdo,istde
*
      OPEN (NFILE,
     :      FILE   = FILNAM,
     :      FORM   = RFORM ,
     :      STATUS = 'UNKNOWN',
     :      IOSTAT = IOS   )
*
      IF (IOS .EQ. 0) THEN
         REWIND (NFILE)
         IERR = 0
      ELSE
         LOC =  LEN_TRIM (FILNAM)
         WRITE (istde,*) 'OPENFL: Error opening file ',FILNAM(1:LOC),
     :            ' as ',RSTAT,';'
         WRITE (istde,*) 'The argument RSTAT=',RSTAT,' is not used !'
         IERR = 1
      ENDIF
*
      RETURN
      END
