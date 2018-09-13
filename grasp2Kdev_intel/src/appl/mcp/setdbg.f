************************************************************************
*                                                                      *
      SUBROUTINE SETDBG (DEBUG, fullname)
*                                                                      *
*   This routine opens the  .dbg  file and sets the arrays that con-   *
*   trol debug printout from the GENMCP program.                       *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, OPENFL.                                *
*                                                                      *
*   Written by Farid A Parpia               Last update: 08 Dec 1992   *
*   Modified by Xinghong He                 Last update: 03 Jul 1998   *
*
*   Used by mcpvu, mcpmpivu
*                                                                      *
************************************************************************
*
      LOGICAL DEBUG,GETYN,LDBPA,LDBPG
      CHARACTER*(*) fullname
      CHARACTER(LEN = LEN (fullname)) FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEFAULT/NDEF

      COMMON/iounit/istdi,istdo,istde
*
*   Initialise the arrays that control the debug printout
*
      DO I = 1, 5
         LDBPA(I) = .FALSE.
         LDBPG(I) = .FALSE.
      ENDDO
*
      IF (NDEF .EQ. 0) THEN
         RETURN
      ENDIF

      WRITE (istde,*) 'Generate debug printout?'
      DEBUG = GETYN ()
      IF (DEBUG) THEN
*
*   The  .dbg  file is formatted; open it on unit 99
*
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         WRITE (istde,*) 'File ', fullname,' will be created as the'
     &,                 ' GENMCP DeBuG Printout File; '
         WRITE (istde,*) ' enter another file name if this is not '
     &,                 'acceptable; null otherwise:'
         READ (*,'(A)') FILNAM
*
         IF (LEN_TRIM (FILNAM) .EQ. 0) FILNAM = fullname
*
    4    CALL OPENFL (99, FILNAM, FORM, STATUS, IERR)
         IF (IERR .NE. 0) THEN
    5       WRITE (istde,*) 'Enter a name for the GENMCP DeBuG Printout'
     &,                    ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LEN_TRIM (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         WRITE (istde,*) ' Print out the machine constants used?'
         LDBPG(1) = GETYN ()
*
*   Set options for angular modules
*
         WRITE (istde,*) ' Printout from LODCSL? (Not used)'
         LDBPA(1) = GETYN ()
         WRITE (istde,*) ' Print out T coefficients?'
         LDBPA(2) = GETYN ()
         WRITE (istde,*) ' Print out Coulomb V coefficients?'
         LDBPA(3) = GETYN ()
         WRITE (istde,*) ' Print out sparse matrix definition arrays?'
         LDBPA(4) = GETYN ()

      ENDIF
*
      RETURN
      END
