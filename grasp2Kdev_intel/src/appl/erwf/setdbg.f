************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
*                                                                      *
*   Written by Farid A Parpia               Last update: 15 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN,LDBPA,LDBPG,LDBPR,YES
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)
      COMMON/DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde
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

      WRITE(istde,*) 'Generate debug printout?'
      YES = GETYN ()
      IF (YES) THEN
*
*   The  .dbg  file is formatted; open it on unit 99
*
         DEFNAM = 'erwf.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         WRITE(istde,*) 'File  erwf.dbg  will be created as the ',
     &                  'ERWF DeBuG Printout File; '
         WRITE(istde,*) 'enter another file name if this is not ',
     &                  'acceptable; null otherwise:'
         READ (*,'(A)') FILNAM
*
         IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       WRITE(istde,*) 'Enter a name for the ERWF DeBuG Printout '
     &                   , 'file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         WRITE(istde,*) ' Print out the machine constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(1) = .TRUE.
         WRITE(istde,*) ' Print out the physical constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(2) = .TRUE.
*
*   Set options for radial modules
*
         WRITE(istde,*) ' Printout from RADGRD?'
         YES = GETYN ()
         IF (YES) LDBPR(1) = .TRUE.
         WRITE(istde,*) ' Printout from NUCPOT?'
         YES = GETYN ()
         IF (YES) LDBPR(2) = .TRUE.
         WRITE(istde,*) ' Printout from TFPOT?'
         YES = GETYN ()
         IF (YES) LDBPR(26) = .TRUE.
*
*   Set options for angular modules
*
         WRITE(istde,*) ' Printout from LODCSL?'
         YES = GETYN ()
         IF (YES) LDBPA(1) = .TRUE.
*
      ENDIF
*
      RETURN
      END
