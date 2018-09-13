************************************************************************
*                                                                      *
      SUBROUTINE SETCSLAjb (NAME,NCORE)
*                                                                      *
*   Open, check, load data from and close the  .csl  file. This file   *
*   is always attached to stream 21.                                   *
*                                                                      *
*   Call(s) to: [RCI92]: LODCSL, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL FOUND
      CHARACTER*24 NAME
      CHARACTER*256 FILNAM
      CHARACTER*15 RECORD
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      COMMON/iounit/istdi,istdo,istde
*
*   The  .csl  file is FORMATTED; it must exist
*
      K=INDEX(NAME,' ')
      FILNAM  = NAME(1:K-1)//'.c'
      FORM = 'FORMATTED'
      STATUS = 'OLD'

*
      CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         WRITE(istde,*) 'Error when opening',FILNAM
         STOP
      ENDIF
*
cbieron setNCF
cbieron Determine NCF from the  .csl  file
*
      CALL setNCF 
*
*   Reposition .csl file
      rewind (21,ERR=911)

*   Check the first record of the file; if not as expected, try again
*
      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.
     :    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         WRITE(istde,*) 'Not a Configuration Symmetry List File;'
         CLOSE (21)
      ENDIF
*
*   Load data from the  .csl  file
*
      CALL LODCSLjb (NCORE)
*
*   Close the  .csl  file
*
      CLOSE (21)
*
      RETURN
 911     WRITE(istde,*) 'Error when rewinding ',FILNAM
         STOP
      END
