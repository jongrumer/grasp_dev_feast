************************************************************************
*                                                                      *
      SUBROUTINE SETISO (fname)
      IMPLICIT INTEGER (I-N)
*                                                                      *
*   Open, check, load data from and close the  .iso  file. This file   *
*   is always attached to stream  22 in  RSCF92,  RCI92,  HFS92, and   *
*   OSCL92.                                                            *
*   Filename - fname is moved to the argument list. This subroutine
*   is ready to replace the one in lib/lib92, but need to check
*   other application programs before doing do.
*                                                                      *
*   Call(s) to: [LIB92]: LODISO, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*   Modified    by Xinghong He            Last revision:  1 Jun 1998   *
*                                                                      *
************************************************************************
      CHARACTER*(*) fname
      CHARACTER(LEN=*), PARAMETER:: FORM = 'FORMATTED', 
     &                              STATUS = 'OLD',
     &                              SIGNATURE = 'Atomic number:',
     &                              MYNAME = 'SETISO'
      CHARACTER(LEN=LEN(SIGNATURE)) str

      LOGICAL FOUND
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
      INQUIRE (FILE = fname, EXIST = FOUND)
      IF (.NOT. FOUND) THEN
         lenf = LEN_TRIM (fname)
         WRITE (istde,*) MYNAME, '- file: ', fname(1:lenf),
     &                           ' does not exist'
         STOP
      ENDIF
*
      CALL OPENFL (22, fname, FORM, STATUS, IERR)
      IF (IERR .NE. 0) THEN
         WRITE (istde,*) 'Error opening isodata file: ', fname(1:lenf) 
         STOP
      ENDIF
*
*   Check the first record of the file; if not as expected, try again
*
      READ (22,'(A)', IOSTAT = IOS) str
      IF ((IOS .NE. 0) .OR.
     :    (str .NE. SIGNATURE)) THEN
         WRITE (istde,*) 'Not an ISOtope Data File;'
         CLOSE (22)
         STOP
      ENDIF
*
*   Load data from the .iso file and then close it.
*
      CALL LODISO
      CLOSE (22)

      RETURN
      END
