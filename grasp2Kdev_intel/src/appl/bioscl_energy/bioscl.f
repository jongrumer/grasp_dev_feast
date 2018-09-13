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
      PRINT *
      !PRINT *, ' Default settings?'
      !YES = GETYN ()
      !PRINT *
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
      PRINT *
      PRINT *, ' Input from a CI calculation?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
        INPCI = 0
      ELSE
        INPCI = 1
      ENDIF
*Rasa -- start
      LDBPR = .false.
      PRINT *, ' Generate debug output?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
        LDBPR(18) = .true.
        PRINT *, ' Enter the cutoff'
        read *, cutoff
      ENDIF
*Rasa -- end
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
      CALL STOPTIME (ncount1, 'BIOSCL3')
CGG      PRINT *
CGG      PRINT *, 'BIOSCL3: Execution complete.'
*
      STOP
      END
