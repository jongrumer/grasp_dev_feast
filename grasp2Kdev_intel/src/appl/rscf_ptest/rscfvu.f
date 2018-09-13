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

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER host*(MPI_MAX_PROCESSOR_NAME), idstring*3

! Things for timing
      INTEGER   ncount1
! cpath uses
      INTEGER       lenperm, lentmp
      CHARACTER*128 startdir, permdir, tmpdir, file_rcsl, file1, file2
!-----------------------------------------------------------------------

*=======================================================================
*  Start mpi --- get processor info: myid, nprocs, host name; and print
*=======================================================================
      startdir = '  '  ;    file_rcsl = '  '
      permdir = '  '   ;    file1     = '  '
      tmpdir = '  '    ;    file2     = '  '
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &                     startdir, permdir, tmpdir, 'RMCDHF_MPI')
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

*=======================================================================

!-----------------------------------------------------------------------
!     myid = 0
!     nprocs = 1

!     write(*,*)
!     write(*,*) 'RMCDHF'
!     write(*,*) 'This program determines the radial orbitals   '
!     write(*,*) 'and the expansion coefficients of the CSFs         '
!     write(*,*) 'in a self-onsistent field proceedure               '
!     write(*,*) 'Input file:  isodata, rcsf.inp, rwfn.inp, mcp.30, ...'
!     write(*,*) 
!    :'Outputfiles: rwfn.out, rmix.out, rmcdhf.sum, rmcdhf.log'
!     write(*,*)

!     CALL starttime (ncount1, 'RMCDHF')


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
      CALL MPI_Bcast (NDEF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)


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

      if(myid == 0) CALL SETSUM ('rmcdhf.sum')

      CALL SETMCP (ncore, nblk0, idblk, 'mcp')
      CALL SETCSL ('rcsf.inp', ncore1, idblk)
            IF (ncore .NE. ncore1) STOP 'rscfvu: ncore'
      CALL MPI_BARRIER (MPI_COMM_WORLD,ierr)

*=======================================================================
*  Gather all remaining information and perform some setup. This
*  part (routine) asks for user-inputs.
*=======================================================================

      CALL GETSCDpmpi (EOL, idblk, 'isodata', 'rwfn.inp')

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

!     CALL stoptime (ncount1, 'RMCDHF')
      CALL stopmpi2 (myid, nprocs, host, lenhost,
     &                     ncount1, 'RMCDHF_MPI')


      STOP
      END
