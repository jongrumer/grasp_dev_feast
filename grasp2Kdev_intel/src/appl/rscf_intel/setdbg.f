************************************************************************
*                                                                      *
      SUBROUTINE SETDBG (dbgfile)
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: OPENFL.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 10 Dec 1992   *
*   Modified bu Xinghong He                 Last update: 06 Jul 1998   *
*                                                                      *
************************************************************************
      CHARACTER*(*) dbgfile
      CHARACTER*(*), PARAMETER:: FORM = 'FORMATTED', status = 'NEW'
      CHARACTER (LEN = LEN (dbgfile)) filnam

      LOGICAL    GETYN, LDBPA, LDBPG, LDBPR
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)
     :      /DEFAULT/NDEF
      COMMON/iounit/istdi,istdo,istde
*
*   Initialise the arrays that control the debug printout
*   These serve as the default settings.
*
      DO I = 1, 5
         LDBPA(I) = .FALSE.
         LDBPG(I) = .FALSE.
      ENDDO

      DO I = 1, 30
         LDBPR(I) = .FALSE.
      ENDDO
*
      IF (NDEF .EQ. 0)  RETURN

*   Even in non-default, the user can choose not to have debug 
*   print-out

      WRITE (istde,'(A)',ADVANCE='NO') 'Generate debug output?  (y/n) '
      IF (.NOT. GETYN ()) RETURN

      lendbg = LEN_TRIM (dbgfile)

      WRITE (istde,*) 'File  ', dbgfile(1:lendbg), '  will be 
     &                created as the RSCF92 DeBuG Printout File;'
      WRITE (istde,*) 'enter another file name if this is not'
     &                , ' acceptable; null otherwise:'

  123 READ (*,'(A)') filnam

      filnam = ADJUSTL (filnam)
      lenfil = LEN_TRIM (filnam)
      IF (lenfil .EQ. 0) THEN
         filnam = dbgfile
      ELSE IF (lenfil .GT. lendbg) THEN
         WRITE (istde,*) 'File name too long, (> ', lendbg, '); redo...'
         GOTO 123
      ENDIF

      CALL OPENFL (99, filnam, FORM, STATUS, IERR)
      IF (IERR .NE. 0) THEN
         WRITE (istde,*) 'File name not accepted; redo...'
         GOTO 123
      ENDIF
*
*   Set options for general printout
*
      WRITE (istde,*) 'Print out the machine constants used?'
      LDBPG(1) = GETYN ()
      WRITE (istde,*) 'Print out the physical constants used?'
      LDBPG(2) = GETYN ()
      WRITE (istde,*) 'Printout from FNDBLK?'
      LDBPG(3) = GETYN ()
      WRITE (istde,*) 'Print out the Hamiltonian matrix?'
      LDBPG(4) = GETYN ()
      WRITE (istde,*) 'Print out the eigenvectors?'
      LDBPG(5) = GETYN (); LDBPG(1:5) = .true.
*
*   Set options for printout from radial modules
*
      WRITE (istde,*) 'Printout from RADGRD?'
      LDBPR(1) = GETYN ()
      WRITE (istde,*) 'Printout from NUCPOT?'
      LDBPR(2) = GETYN ()
      WRITE (istde,*) 'Printout from LODRWF?'
      LDBPR(3) = GETYN ()
      WRITE (istde,*) 'Print out I(ab) integrals?'
      LDBPR(4) = GETYN ()
      WRITE (istde,*) 'Print out Slater integrals?'
      LDBPR(10) = GETYN ()
      WRITE (istde,*) 'Make summary printout on progress'
     &, ' of each iteration in SOLVE?'
      LDBPR(22) = GETYN ()
      WRITE (istde,*) 'Tabulate and make printer plots'
     &, ' of subshell radial functions on'
     &, ' each iteration in SOLVE?'
      LDBPR(23) = GETYN ()
      WRITE (istde,*) 'Tabulate and make printer plots'
     &, ' of subshell radial functions'
     &, ' after each SCF cycle?'
      LDBPR(24) = GETYN ()
      WRITE (istde,*) 'Tabulate and make printer plots'
     &, ' of subshell radial functions on'
     &, ' convergence?'
      LDBPR(25) = GETYN ()
      WRITE (istde,*) 'List compositions of exchange'
     &, ' potentials?'
      LDBPR(27) = GETYN ()
      WRITE (istde,*) 'Tabulate and make printer plots'
     &, ' of exchange potentials?'
      LDBPR(28) = GETYN ()
      WRITE (istde,*) 'List compositions of direct'
     &, ' potentials?'
      LDBPR(29) = GETYN ()
      WRITE (istde,*) 'Tabulate and make printer plots'
     &, ' of direct potentials?'
      LDBPR(30) = GETYN (); LDBPR(1:30) = .true. 
*
*   Set options for printout of angular coefficients
*
      WRITE (istde,*) ' Printout from LODCSL?'
      LDBPA(1) = GETYN ()
      WRITE (istde,*) ' Print out T coefficients?'
      LDBPA(2) = GETYN ()
      WRITE (istde,*) ' Print out V coefficients?'
      LDBPA(3) = GETYN (); LDBPA(1:3) = .true.

      RETURN
      END
