************************************************************************
*
      SUBROUTINE SETCSH (nfile, name, ncore)
*
*   Open, check the CSL file and load the load (via lodcsh) data from 
*   the header lines. It is designed to replace all kinds of "setcsl" 
*   routines within GRASP packages.
*
*   Routines called: lodcsh
*
*   ncore is output parameter
*
*   Written by Xinghone He                Last revision: 23 Dec 1997 
*
************************************************************************
      CHARACTER name*(*)
      COMMON/iounit/istdi,istdo,istde

*     ...Locals
      LOGICAL FOUND
      CHARACTER*15 RECORD

      length = LEN_TRIM(name)

      INQUIRE (FILE = name, EXIST = FOUND)
      IF (.NOT. FOUND) THEN
         WRITE (istde,*) name(1:length),' does not exist'
         STOP
      ENDIF

      INQUIRE (UNIT = nfile, OPENED = FOUND)
      IF (FOUND) THEN
         WRITE (istde,*) 'Unit ', nfile, ' has been used elsewhere'
         STOP
      ENDIF

      OPEN (nfile, FILE = name, STATUS = 'OLD', IOSTAT = IOS)
      IF (IOS .NE. 0) THEN
         WRITE (istde,*) 'Error when opening ', name
         STOP
      ENDIF
*
*   Check the first record of the file; if not as expected, try again
*
      READ (nfile, '(1A15)', IOSTAT = IOS) RECORD

      IF (IOS .NE. 0 .OR. RECORD(1:15) .NE. 'Core subshells:') THEN
         WRITE (istde,*) 'Not a Configuration Symmetry List File;'
         CLOSE (nfile)
         STOP
      ENDIF
*
*   Load data from the  .csl  file
*
      CALL LODCSH (nfile, NCORE)

      RETURN
      END
