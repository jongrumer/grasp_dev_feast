************************************************************************
*                                                                      *
      PROGRAM BIOSCL
*                                                                      *
*   This program calculates the transition parameters for a            *
*   transition between an initial and a final state                    *
*   The program assumes that the radial orbitals of the initial        *
*   and final state have been transformed by the BIOTRA program        *
*   as to become biorthonormal, in which case the normal Racah         *
*   algebra can be used.                                               *
*                                                                      *
*   Written by Per Jonsson,   Department of Computer Science           *
*                             Vanderbilt University, USA               *
*                                                                      *
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN,YES,LDBPR
      CHARACTER*24 NAME(2)
      DOUBLE PRECISION CUTOFF

      COMMON/DEFAULT/NDEF,NDUMP,INPCI
      COMMON/DEBUGR/LDBPR(30), CUTOFF


      NTEST = 1001
      CALL STARTTIME (ncount1, 'BIOSCL3')
CGG      PRINT *
CGG      PRINT *, ' BIOSCL3: Execution begins ...'
      write(*,*)
      write(*,*) 'RTRANSITION'
      write(*,*) 'This program computes transition parameters from'
      write(*,*) 'transformed wave functions'
      write(*,*) 'Input files:  isodata, name1.c, name1.bw, name1.(c)bm'
      write(*,*) '              name2.c, name2.bw, name2.(c)bm         '
      write(*,*) '              optional, name1.lsj.lbl, name2.lsj.lbl'
      write(*,*) '              name1.name2.KT (optional angular files)'
      write(*,*) 'Output files: name1.name2.(c)t                       '
      write(*,*) '              optional, name1.name2.(c)t.lsj         '
      write(*,*) '              name1.name2.KT (angular files)         '
      write(*,*) 'Here K is parity and rank of transition: -1,+1 etc   '

      PRINT *
      PRINT *, ' Default settings?'
      YES = GETYN ()
      YES = .TRUE.
      IF (YES) THEN
        NDEF = 0
        NDUMP = 1
      ELSE
        NDEF = 1
        PRINT *, ' Dump angular data to file?'
        YES = GETYN ()
        IF (YES) THEN
          NDUMP = 1
        ELSE
          NDUMP = 0
        ENDIF
      ENDIF
*Rasa -- start
      LDBPR = .false.
C      PRINT *, ' Generate debug output?'
C      YES = GETYN ()
C      IF (YES) THEN
C        LDBPR(18) = .true.
C        PRINT *, ' Enter the cutoff'
C        read *, cutoff
C      ENDIF
*Rasa -- end
      PRINT *, ' Input from a CI calculation?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
        INPCI = 0
      ELSE
        INPCI = 1
      ENDIF
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Obtain the names of the initial and final state files
*
      CALL FNAME(NAME)
*
*   Open, check, load data from, and close, the initial and final state
*   .csl  files. These files are then merged to one file.
*
      CALL MRGCSL(NAME)
*
*   Open, check, load data from, and close, the merged .csl  file
*
      CALL SETCSLM
*
*   Read mixing coefficients
*
C      CALL READMIX(NAME,INPCI)
*
*   Test mixing coefficients
*
C      IF (NTEST.GT.1000) CALL TESTMIX
*
*   Get the remaining information
*
      CALL GETOSD(NAME)
*
*   Open and append a summary of the inputs to the  .sum  file
*
      ILBL = 0
      CALL STRSUM(NAME,INPCI,ILBL)
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the transition calculation
*
      CALL OSCL(NAME)
*
*   Print completion message
*
      CALL STOPTIME (ncount1, 'RTRANSITION')
CGG      PRINT *
CGG      PRINT *, 'BIOSCL3: Execution complete.'
*
      STOP
      END
