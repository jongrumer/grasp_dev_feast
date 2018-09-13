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
      CHARACTER*24 NAME
      CHARACTER*256 FILNAM
      CHARACTER*15 RECORD
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
*   The  .csl  file is FORMATTED; it must exist
*
      K=INDEX(NAME,' ')
c     K=LEN_TRIM(NAME)+1
      FILNAM  = TRIM(NAME)//'.c'
      FORM = 'FORMATTED'
      STATUS = 'OLD'

*
c     CALL OPENFL (21,TRIM(FILNAM),FORM,STATUS,IERR)
c     IF (IERR .EQ. 1) THEN
c        PRINT *, 'Error when opening',FILNAM
c        STOP
c     ENDIF
C      print *,filnam,'filnam',name
      open(21,file=TRIM(NAME)//'.c',form='formatted',status='old')
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
      CALL LODCSLO (NCORE)
*
*   Close the  .csl  file
*
      CLOSE (21)
*
      RETURN
      END
