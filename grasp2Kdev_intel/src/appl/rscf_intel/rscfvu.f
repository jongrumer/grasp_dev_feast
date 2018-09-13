************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***       ******    *****    *****   *******   *****    *****        ***
***       **   **  **   **  **   **  **       **   **  **   **       ***
***       **   **  **       **       **       **   **       **       ***
***       ******    *****   **       ****      *****       **        ***
***       **  **        **  **       **          **       **         ***
***       **   **  **   **  **   **  **         **      **           ***
***       **   **   *****    *****   **        **      *******       ***
***                                                                  ***
***            Relativistic Self-Consistent-Field Program            ***
***                                                                  ***
***   This program is a derivative of GRASP2 (F. A. Parpia, I. P.    ***
***   Grant, and C. F. Fischer, 1990)                                ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***   Modified by C.F. FIscher for block input                       ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM RSCFVU
      IMPLICIT NONE
*                                                                      *
*   Entry routine for RSCFVU. Controls the entire computation.      *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RSCF92]: CHKPLT, setcsl, setdbg
*                        SETMIX, SETRES, SETSUM, STRSUM.
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
*   Block version by Xinghong He          Last revision: 17 Aug 1998   *
*                                                                      *
************************************************************************
*
CGG      INTEGER, PARAMETER::nblk0 = 20	! Maximum number of blocks
      INTEGER, PARAMETER::nblk0 = 50	! Maximum number of blocks
      CHARACTER*7 idblk(nblk0)

      LOGICAL EOL, GETYN, YES

      INTEGER NDEF, NCORE, istdi, istdo, istde, ncore1
      COMMON/DEFAULT/NDEF
      COMMON/CORE/NCORE
      COMMON/iounit/istdi,istdo,istde

      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr

! Things for timing
      INTEGER   ncount1
!-----------------------------------------------------------------------
      myid = 0
      nprocs = 1

      write(*,*)
      write(*,*) 'RMCDHF'
      write(*,*) 'This program determines the radial orbitals   '
      write(*,*) 'and the expansion coefficients of the CSFs         '
      write(*,*) 'in a self-onsistent field proceedure               '
      write(*,*) 'Input file:  isodata, rcsf.inp, rwfn.inp, mcp.30, ...'
      write(*,*) 
     :'Outputfiles: rwfn.out, rmix.out, rmcdhf.sum, rmcdhf.log'
      write(*,*)

      CALL starttime (ncount1, 'RMCDHF')


      OPEN(UNIT=734,FILE='rmcdhf.log',STATUS='UNKNOWN')

*=======================================================================
*  Get NDEF
*=======================================================================

      IF (myid .EQ. 0) THEN
         WRITE (istde,'(A)',ADVANCE='NO') ' Default settings?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
            WRITE(734,'(A)') 'y            ! Default settings'
         ELSE
            NDEF = 1
         ENDIF
      ENDIF

*=======================================================================
*
*  Checks and settings... Mostly done in backyard.
*
*    CHKPLT - compatibility of plant substitutions
*    SETDBG - Debug output control parameters
*    SETMC - machine- and installation- dependent constants
*    SETCON - physical constants
*    SETSUM - open the summary file
*    SETMCP - open and check the  .mcp  files
*    SETCSL - open, check, load data from, and close the  .csl  file
*    STRSUM - append a summary of the inputs to the  .sum  file
*    SETMIX - mixing coefficients file setup
*    FACTT - table of logarithms of factorials setup
*=======================================================================

      CALL CHKPLT ('RSCFVU')
      CALL SETDBG ('rscf92.dbg')
      CALL SETMC
      CALL SETCON

      CALL SETSUM ('rmcdhf.sum')

      CALL SETMCP (ncore, nblk0, idblk, 'mcp')
      CALL SETCSL ('rcsf.inp', ncore1, idblk)
            IF (ncore .NE. ncore1) STOP 'rscfvu: ncore'

*=======================================================================
*  Gather all remaining information and perform some setup. This
*  part (routine) asks for user-inputs.
*=======================================================================

      CALL GETSCD (EOL, idblk, 'isodata', 'rwfn.inp')

      IF (myid .EQ. 0) THEN
         CALL STRSUM
         IF (EOL) CALL SETMIX ('rmix.out')
      ENDIF

      CALL FACTT

*=======================================================================
*  Proceed with the SCF calculation close all files except
*  the  .sum  file
*=======================================================================

      CALL scf (EOL, 'rwfn.out')

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================

      CALL stoptime (ncount1, 'RMCDHF')

      STOP
      END
