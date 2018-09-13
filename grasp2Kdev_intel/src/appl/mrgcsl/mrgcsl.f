************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***      **    **  ******    *****    *****    *****   **            ***
***      ***  ***  **   **  **   **  **   **  **   **  **            ***
***      ** ** **  **   **  **       **       **       **            ***
***      ** ** **  ******   **  ***  **        *****   **            ***
***      **    **  **  **   **   **  **            **  **            ***
***      **    **  **   **  **   **  **   **  **   **  **            ***
***      **    **  **   **   *****    *****    *****   *******       ***
***                                                                  ***
***        Program to merge Configuration Symmetry List Files        ***
***                                                                  ***
***                         GRASP92 Version                          ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM MRGCSL
*                                                                      *
*   Entry routine for MRGCSL. Controls the entire computation.         *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC.                                        *
*               [MRGCSL]: CHKPLT, LDCSL1, LDCSL2, MERG12, SETDBG.      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      WRITE(*,*) 
      WRITE(*,*) 'RCSFMERGE'
      WRITE(*,*) 'Merges two files of CSFs and weeds out duplicates'
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
*   Merge the two  .csl  lists, eliminating any repeated CSFs
*
      CALL MERG12 (NCORER,NCORE)
*
*   Split into block
*
      CALL RCSFBLOCK

*
*   Print completion message
*
      


      PRINT *, 'RCSFMERGE: Execution complete.'
*
      STOP
      END
