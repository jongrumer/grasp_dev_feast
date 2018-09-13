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
*                                                                      *
************************************************************************
*
      LOGICAL GETYN,YES
      COMMON/DEFAULT/NDEF
      EXTERNAL CONSTS
      COMMON/CONS/ZERO,HALF,TENTH,ONE,TWO,THREE,TEN
      COMMON/iounit/istdi,istdo,istde
      CHARACTER*1 answer
*   Startup message
*
      WRITE (istde,*) 'ERWF: Execution begins ...'
      WRITE (istde,*) 'Estimating Relativistic Wave Functions: ',
     &                ' Output file = rwfn.inp'
*
      WRITE (istde,*) 'Default settings ?'
      YES = GETYN ()
      IF (YES) THEN
         NDEF = 0
      ELSE
         NDEF = 1
      ENDIF

      WRITE (istde,*) 'Normal or reversed orbital order (n/r)?'
      READ(*,'(a)') answer
      IF ((answer.eq.'n').or.(answer.eq.'N')) THEN
         NORDER = 1
      ELSE IF ((answer.eq.'r').or.(answer.eq.'R')) THEN
         NORDER = 0
      ELSE
         WRITE(*,*) 'Expected answer n or r'
         STOP
      END IF
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
      CALL SETCSH (21, 'rcsl.inp', ncore)

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
*   norder = 1 normal orbital order, orthogonalize against orbitals in
*              the beginning of the orbital list
*   norder = 0 reversed orbital order, orthogonalize against orbitals 
*              at the end of the orbital list
      IF (NORDER.EQ.1) THEN
         CALL ORTHSC
      ELSE
         CALL ORTHSCR
      END IF
*
*   Write the subshell radial wavefunctions out
*
      CALL WRTRWF
*
*   Print completion message
*
      WRITE (istde,*) 'ERWF: Execution complete.'
*
      STOP
      END
