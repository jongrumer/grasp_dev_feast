************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***               *******  *******   **    **  *******               ***
***               **       **    **  **    **  **                    ***
***               **       **    **  **    **  **                    ***
***               *****    *******   ** ** **  *****                 ***
***               **       **  **    ** ** **  **                    ***
***               **       **   **   ** ** **  **                    ***
***               *******  **    **   **  **   **                    ***
***                                                                  ***
***   Program to create a  .rwf  file for RSCF92 by merging one or   ***
***   more existing  .rwf  files or generating Thomas-Fermi or hy-   ***
***   drogenic functions.                                            ***
***                                                                  ***
***                              GRASP92                             ***
***           F. A. Parpia, C. F. Fischer, and I. P. Grant           ***
***                                                                  ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM ERWF
      IMPLICIT REAL*8          (A-H, O-Z)
*                                                                      *
*   Entry routine for RCI92. Controls the entire computation.          *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RCI92]: CHKPLT, MATRIX, SETCSL, SETDBG, SETMIX,       *
*                        SETRES, SETSUM, STRSUM.                       *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 15 Dec 1992   *
*   Modified by Per Jonsson to always have orbitals in                 *
*   canonical order             March 2014                             *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN,YES
      COMMON/DEFAULT/NDEF
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde

*   Startup message
*
      write(*,*) 
      write(*,*) 'RWFNESTIMATE'
      write(*,*) 'This program estimates radial wave functions '
      write(*,*) 'for orbitals'
      write(*,*) 'Input files: isodata, rcsf.inp, optional rwfn file'
      write(*,*) 'Output file: rwfn.inp'
      write(*,*)
*

      WRITE (istde,*) 'Default settings ?'
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
      IF (NDEF.NE.0) CALL SETSUM
*
*   Open, check, load data from, and close the  .csl  file
*
      CALL SETCSH (21, 'rcsf.inp', ncore)

*   Arrange orbitals in canonical order
 
      CALL CANONICAL
*
*   Hydrogenic screen parameters for all orbitals
*
      CALL screenpar (ncore)
*
*   Determine other relevant information
*
      CALL GETINF
*
*   Write the first part of the  .sum  file
*
      IF (NDEF.NE.0) CALL STRSUM
*
*   Generate the subshell radial wavefunctions
*
      CALL GENRWF
*
*   Orthogonalize the radial orbitals
*
      CALL ORTHSC
*
*   Write the subshell radial wavefunctions out
*
      CALL WRTRWF
*
*   Print completion message
*
      WRITE(*,*)
      WRITE (istde,*) 'RWFNESTIMATE: Execution complete.'
*
      STOP
      END
