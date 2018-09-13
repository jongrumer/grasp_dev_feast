************************************************************************
*                                                                      *
      PROGRAM JJREDUCE
*                                                                      *
*   From a set of CSLs this program identifies the ones that           *
*   interact with a given multireference                               *
*                                                                      *
*   This program is a slight modification of the GENMCP program        *
*                                                                      *
*   Written by Per Jonsson                                             *
*                                                                      *
*   Modified by Gediminas Gaigalas for new spin-angular integration.   *
*                                                                      *
************************************************************************
*
      LOGICAL RESTRT,DEBUG,GETYN,YES
      COMMON/DEFAULT/NDEF
*
      PRINT *
      PRINT *, 'JJREDUCE: Execution begins ...'
      PRINT *

      PRINT *, 'Default settings?'
      YES = GETYN ()
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF


*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG (DEBUG)
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Open the  .sum  file
*
C      IF (NDEF.NE.0) CALL SETSUM
*
*   Cetermine the names of the file to be reduce and the mr file
*
*   Open, check, load data from, and close the  csl  file
*
      CALL SETCSLB
*
*   Determine where the multireference CSFs reside in the csl file
*
      CALL IDENTY
*
*   Set up the  mcp  files; determine if this is a restart
*
C      CALL SETMCP(RESTRT)
      RESTRT = .FALSE.
*
*   Append a summary of the inputs to the  .sum  file
*
C      IF (NDEF.NE.0) CALL STRSUM
*
*   Set up the table of logarithms of factorials for use by
*   angular modules
*
      CALL FACTT
*
*   Proceed with the generation of MCP coefficients
*
      CALL MCP(RESTRT)
*
*   All done; close files that are still open
*
      CLOSE (24)
C      IF (DEBUG) CLOSE (99)
*
*   Print completion message
*
      PRINT *
      PRINT *, 'JJREDUCE: Execution complete.'
      PRINT *
      PRINT *, ' The reduced list is in rcsl.out'
*
      STOP
      END
