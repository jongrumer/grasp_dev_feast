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
      PROGRAM RSCFMPIVU
      IMPLICIT NONE
*                                                                      *
*   Entry routine for RSCFMPIVU. Controls the entire computation.      *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC, SETCON.                                *
*               [RSCF92]: CHKPLT, setcslmpi, setdbgwrap, 
*                        SETMIX, SETRES, SETSUM, STRSUM.
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 31 Dec 1992   *
*   MPI version by Xinghong He            Last revision: 05 Aug 1998   *
*                                                                      *
************************************************************************
*
      INTEGER, PARAMETER::nblk0 = 20	! Maximum number of blocks
      CHARACTER*8 idblk(nblk0)

      LOGICAL EOL, GETYN, YES

      INTEGER NDEF, NCORE, istdi, istdo, istde, ncore1
      COMMON/DEFAULT/NDEF
      COMMON/CORE/NCORE
      COMMON/iounit/istdi,istdo,istde

! mpi common is written here and kept unchanged through out.
      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr, lenhost
      COMMON /mpi/ myid, nprocs, ierr
      CHARACTER host*(MPI_MAX_PROCESSOR_NAME), idstring*3

! Things for timing
      INTEGER   ncount1

! cpath uses
      INTEGER       lenperm, lentmp
      CHARACTER*128 startdir,permdir,tmpdir,file_rcsl,file1,file2
!-----------------------------------------------------------------------

*=======================================================================
*  Start mpi --- get processor info: myid, nprocs, host name; and print
*=======================================================================
      startdir = '  '  ;    file_rcsl = '  '
      permdir = '  '   ;    file1     = '  '
      tmpdir = '  '    ;    file2     = '  '
      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &    startdir, permdir, tmpdir, 'RSCFMPIVU',lentmp,lenperm)
      WRITE (idstring, '(I3.3)') myid
      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

*=======================================================================
*  Get NDEF on node-0 and then send to all nodes   
*=======================================================================

      IF (myid .EQ. 0) THEN
         WRITE (istde,*)
         WRITE (istde,'(A)',ADVANCE='YES') 'Default settings?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
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
*    SETDBGmpi - Debug output control parameters
*    SETMC - machine- and installation- dependent constants
*    SETCON - physical constants
*    SETSUM - open the log file
*    SETMCP - open and check the  .mcp  files
*    SETCSLmpi - open, check, load data from, and close the  .csl  file
*    STRSUM - append a summary of the inputs to the  .logfile
*    SETMIX - mixing coefficients file setup
*    FACTT - table of logarithms of factorials setup
*=======================================================================

      CALL CHKPLT ('RSCFMPIVU')
      CALL SETDBGmpi (permdir(1:lenperm) // '/rscf92.dbg')
      CALL SETMC
      CALL SETCON

      IF (myid .EQ. 0) CALL SETSUM (permdir(1:lenperm-1) // '/rscf.log')

      CALL SETMCP (ncore, nblk0, idblk, 'mcp' // idstring)
      if(myid == 0) file_rcsl = permdir(1:lenperm-1)//'/rcsf.inp'
      CALL SETCSLmpi (file_rcsl, ncore1, idblk)
            IF (ncore .NE. ncore1) STOP 'rscfmpivu: ncore'

*=======================================================================
*  Gather all remaining information and perform some setup. This
*  part (routine) asks for user-inputs.
*=======================================================================
      if(myid == 0) file1 = permdir(1:lenperm-1) // '/isodata'
      if(myid == 0) file2 = permdir(1:lenperm-1) // '/rwfn.inp'
      CALL GETSCDmpi (EOL, idblk,file1,file2) 
!     &            permdir(1:lenperm-1) // '/isodata',
!     &            permdir(1:lenperm-1) // '/rwfn.inp')
      file1 = '  '
      if(myid == 0) file1 = permdir(1:lenperm-1) // '/rmix.out'
      IF (myid .EQ. 0) THEN
         CALL STRSUM
         IF (EOL) CALL SETMIX (file1)
      ENDIF

      CALL FACTT

*=======================================================================
*  Proceed with the SCF calculation close all files except
*  the  .log  file
*=======================================================================
      file1 = '  '
      if(myid == 0) file1 = permdir(1:lenperm-1) // '/rwfn.out'
      CALL scfmpi (EOL, file1)

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost, 
     &                     ncount1, 'RSCFMPIVU')

      STOP
      END
