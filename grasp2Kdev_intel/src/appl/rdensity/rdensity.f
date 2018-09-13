************************************************************************
*                                                                      *
      PROGRAM RDENSITY
*                                                                      *
*   Entry routine for RDENSITY. Controls the entire computation.           *
*                                                                      *
*   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
*               [SMS92]: CHKPLT, GETSMD, SETDBG, SETSUM                *
*                        STRSUM, SETDENS.                                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Per Jonsson                Last revision: 17 Jan 1996   *
*   Modify  by Gediminas Gaigalas                        26 Oct 2009   *
*   Modified by J. Ekman                                 25 Nov 2013   *
*                                                                      *
************************************************************************
*
      DOUBLE PRECISION DR2
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      COMMON/DEFAULT/NDEF

      PRINT *
      PRINT *, 'RDENSITY: Execution begins ...'

      PRINT *
      PRINT *, 'Default settings?'
      YES = GETYN ()
      PRINT *
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

    9 PRINT *, 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         PRINT *, 'Names may not start with a blank'
         GOTO 9
      ENDIF
      PRINT *

      PRINT *, 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
      PRINT *

      ! JE ADD
   10 PRINT *, 'Three- or four-parameter electron density fit within nuc
     :lear volume?'
      PRINT *, '3: b0 + b2*r^2 + b4*r^4'
      PRINT *, '4: b0 + b2*r^2 + b4*r^4 + b6*r^6'
      READ(*,*) NOPAR
      PRINT *
      IF ((NOPAR.NE.3).AND.(NOPAR.NE.4)) THEN
         PRINT *, 'Please answer 3 or 4!'
         GOTO 10
      ENDIF

      ! JE ADD
   11 PRINT *, 'Enter difference in <r^2> value between isotopes (fm^2)'
      READ(*,*) DR2
      PRINT *

*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Set up the physical constants
*
      CALL SETCON
*
*   Open the  .i  file
*
      CALL SETSUM(NAME,NCI)
*
*   Open the  .d  file                                       ! JE ADD
*
      CALL SETDENS(NAME,NCI)                                 ! JE ADD
*
*   Open, check, load data from, and close, the  .csl  file
*
      CALL SETCSLA(NAME,ncore_not_used)
*
*   Get the remaining information
*
      CALL GETSMD(NAME)
*
*   Get the eigenvectors
*
*      PRINT *, 'Block format?'
*      YES = GETYN ()
*      PRINT *

*      IF (YES) THEN
         CALL GETMIXBLOCK(NAME,NCI)
*      ELSE
*         IF (NCI.EQ.0) THEN
*            CALL GETMIXC(NAME)
*         ELSE
*            CALL GETMIXA(NAME)
*         ENDIF
*      ENDIF
*
*   Append a summary of the inputs to the  .sum  file
*
      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the RDENSITY calculation
*
      CALL RDENSITY_CAL(NAME,NOPAR,DR2)
*
*   Print completion message
*
      PRINT *
      PRINT *, 'RDENSITY: Execution complete.'
*
*      STOP
      END
