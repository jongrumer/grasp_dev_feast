************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***                **   **   *****    *****   **                     ***
***                 ** **   **   **  **   **  **                     ***
***                  ***    **       **       **                     ***
***                  ***    **        *****   **                     ***
***                  ***    **            **  **                     ***
***                 ** **   **   **  **   **  **                     ***
***                **   **   *****    *****   *******                ***
***                                                                  ***
***   Utility to delete a set of CSFs from the CSL File              ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM XCSL
*                                                                      *
*   Entry routine for XCSL. Controls the entire computation.           *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC.                                        *
*               [XCSL]: CHKPLT, LDCSL1, LDCSL2, XCLD12, SETDBG.        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 25 Sep 1993   *
*                                                                      *
************************************************************************
*
      WRITE(*,*)
      WRITE(*,*) 'RCSFDELETE'
      WRITE(*,*) 'This program subtracts CSFs in file2 from the CSFs in'
      WRITE(*,*) 'file1. '
      WRITE(*,*) 'Input files: file1 (full name of first CSF file)'
      WRITE(*,*) '             file2 (full name of second CSF file)' 
      WRITE(*,*) 'Output file: rcsf.out'
      WRITE(*,*)
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
*   Load the first  .csl  file
*
      CALL LDCSL1 (NCORER)
*
*   Load the second  .csl  file
*
      CALL LDCSL2 (NCORE)
*
*   Eliminate CSFs common to both lists
*
      CALL XCLD12 (NCORER,NCORE)
*
*   Split into symmetry blocks
*
      CALL RCSFBLOCK

*
*   Print completion message
*
      PRINT *, 'RCSFDELETE: Execution complete.'
*
      STOP
      END
