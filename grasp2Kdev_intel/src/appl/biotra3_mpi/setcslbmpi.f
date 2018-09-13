************************************************************************
*                                                                      *
      SUBROUTINE SETCSLB(NAME,NCORE)
*                                                                      *
*   Open, check, load data from and close the  .csl  file. This file   *
*   is always attached to stream 21.                                   *
*                                                                      *
*   Call(s) to: [RCI92]: LENGTH, LODCSL, OPENFL.                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL FOUND
      CHARACTER*128 NAME
      CHARACTER*256 FILNAM
      CHARACTER*15 RECORD
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
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
         PRINT *, 'Error when opening',FILNAM
         STOP
      ENDIF
*
*   Check the first record of the file; if not as expected, try again
*
      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.
     :    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         PRINT *, 'Not a Configuration Symmetry List File;'
         CLOSE (21)
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
