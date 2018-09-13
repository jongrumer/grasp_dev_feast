************************************************************************
*                                                                      *
      SUBROUTINE SETCSL(NCORE)
*                                                                      *
*   Open, check, load data from and close the  .csl  file. This file   *
*   is always attached to stream 21.                                   *
*                                                                      *
*   Call(s) to: [RCI92]:  LODCSL, OPENFL.                              *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL FOUND
      CHARACTER*256 FILNAM
      CHARACTER*15 RECORD
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      COMMON/iounit/istdi,istdo,istde
*
*   The  .csl  file is FORMATTED; it must exist
*
      DEFNAM = 'rcsf.inp'
      FORM = 'FORMATTED'
      STATUS = 'OLD'
*
*   Look for  csl file
*
      INQUIRE (FILE = DEFNAM, EXIST = FOUND)
*
      IF (FOUND) THEN
         FILNAM = DEFNAM
      ELSE
         WRITE(istde,*) 'rcsf.inp does not exist'
         STOP
      ENDIF
*
*   Open csl file and check the file header
*
      CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE(istde,*) 'Error when opening rcsf.inp'
         STOP
      ENDIF

      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.
     :    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         WRITE(istde,*) 'Not a Configuration Symmetry List File;'
         CLOSE (21)
         STOP
      ENDIF
*
*   Load data from the  .csl  file
*
      CALL LODCSL (NCORE)
*
*   Close the  .csl  file
*
      CLOSE (21)
*
      RETURN
      END
