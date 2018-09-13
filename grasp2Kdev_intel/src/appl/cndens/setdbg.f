************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: OPENFL.                                       *
*                                                                      *
*   Written by Farid A Parpia               Last update: 23 Dec 1992   *
*                                                                      *
************************************************************************
! LENGTH --> LEN_TRIM
! 1997.02.10
      LOGICAL LDBPA,LDBPG,LDBPR
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
      CHARACTER*1 ANSWER
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)
     :      /DEFAULT/NDEF
*
*   Initialise the arrays that control the debug printout
*
      DO 1 I = 1,5
         LDBPA(I) = .FALSE.
    1 CONTINUE
*
      DO 2 I = 1,5
         LDBPG(I) = .FALSE.
    2 CONTINUE
*
      DO 3 I = 1,30
         LDBPR(I) = .FALSE.
    3 CONTINUE
*
      IF (NDEF.EQ.0) THEN
         RETURN
      ENDIF

      PRINT *
      PRINT *, 'Generate debug printout?'
      READ (*,'(A)') ANSWER
      IF (ANSWER .EQ. 'y') THEN
*
*   The  .dbg  file is formatted; open it on unit 99
*
         DEFNAM = 'rci92.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         PRINT *, 'File  rci92.dbg  will be created as the'
         PRINT *, ' RCI92 DeBuG Printout File; enter another'
         PRINT *, ' file name if this is not acceptable;'
         PRINT *, ' <cr> otherwise:'
         READ (*,'(A)') FILNAM
*
         IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       PRINT *, 'Enter a name for the RCI92 DeBuG Printout'
            PRINT *, ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         PRINT *
         PRINT *, ' Print out the machine constants used?'
         READ (*,'(A)') ANSWER
         IF (ANSWER .EQ. 'y') LDBPG(1) = .TRUE.
*
*   Set options for angular modules
*
         PRINT *
         PRINT *, ' Printout from LODCSL?'
         READ (*,'(A)') ANSWER
         IF (ANSWER .EQ. 'y') LDBPA(1) = .TRUE.
*
      ENDIF
*
      RETURN
      END
