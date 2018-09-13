************************************************************************
*                                                                      *
      PROGRAM RIS
*                                                                      *
*   Entry routine for RIS3. Controls the entire computation.           *
*                                                                      *
*   Call(s) to: [LIB92]: GETMIX, SETCSL, SETMC, SETCON.                *
*               [SMS92]: CHKPLT, GETSMD, SETDBG, SETSUM                *
*                        STRSUM.                                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Per Jonsson                Last revision: 17 Jan 1996   *
*   Modify  by Gediminas Gaigalas                        26 Oct 2009   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN, YES
      CHARACTER*24 NAME
      COMMON/DEFAULT/NDEF

      write(*,*)
      write(*,*) 'RIS'
      write(*,*) 'This the isotope shift program '
      write(*,*) 'Input files:  isodata, name.c, name.(c)m, name.w'
      print*,'              name.IOB, name.IOT (optional angular files)'
      write(*,*) 'Output files: name.(c)i '
      print*,'              name.IOB, name.IOT (optional angular files)'
      PRINT *, 'Default settings?'
      YES = GETYN ()
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

   10 PRINT *, 'Name of state'
      READ(*,'(A)') NAME
      K=INDEX(NAME,' ')
      IF (K.EQ.1) THEN
         PRINT *, 'Names may not start with a blank'
         GOTO 10
      ENDIF
      PRINT *, 'Mixing coefficients from a CI calc.?'
      YES = GETYN ()
      IF (YES) THEN
         NCI = 0
      ELSE
         NCI = 1
      ENDIF
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
*   Open the  .sum  file
*
      CALL SETSUM(NAME,NCI)
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
C      CALL STRSUM
*
*   Set up the table of logarithms of factorials
*
      CALL FACTT
*
*   Proceed with the RIS3 calculation
*
      CALL RIS_CAL(NAME)
*
*      STOP
      END
