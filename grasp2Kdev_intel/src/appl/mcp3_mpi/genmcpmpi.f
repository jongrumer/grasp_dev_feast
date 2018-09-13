
************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***        *****   ******  **   **  **    **   *****   ******        ***
***       **   **  **      ***  **  ***  ***  **   **  **   **       ***
***       **       **      ***  **  ** ** **  **       **   **       ***
***       **  ***  ****    ** ****  ** ** **  **       ******        ***
***       **   **  **      **  ***  **    **  **       **            ***
***       **   **  **      **   **  **    **  **   **  **            ***
***        *****   ******  **   **  **    **   *****   **            ***
***                                                                  ***
***   Program for generating the energy expression for H(DC). This   ***
***   program is a derivative of GRASP2 (F. A. Parpia, I. P. Grant,  ***
***   and C. F. Fischer, 1990).                                      ***
***                                                                  ***
***                         GRASP92 Version                          ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
*                                                                      *
*   Entry routine for GENMCP. Controls the entire computation.         *
*                                                                      *
*   Call(s) to: [LIB92]:  GETYN, SETMC, LODCSH2                        *
*               [GENMCP]: CHKPLT, SETDBG, SETSUM, SETCSL, SETMCP,      *
*                         STRSUM, SETTMP, MCP.                         *
*               [NJGRAF]: FACTT.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 11 Dec 1992   *
*   MPI version by Xinghong He            Last revision: 29 Jun 1998   *
*                                                                      *
************************************************************************

      include 'parameters.def'
CGG      PARAMETER (NNNWP = 30)
      PARAMETER (nblk0 = 20)
      CHARACTER*8 idblk(nblk0)

      LOGICAL DEBUG, RESTRT, GETYN, YES

      COMMON/DEFAULT/NDEF

      POINTER (pncfblk, ncfblk(*))
      COMMON/hblock/nblock, pncfblk

      INTEGER*4 IQA,JQSA,JCUPA

      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /STAT/PNTJQS,PNJCUP

      COMMON/MCPA/KMAX
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
      CHARACTER*128 startdir, permdir, tmpdir
!-----------------------------------------------------------------------

      IF (myid .EQ. 0) THEN
        WRITE(*,*)
        WRITE(*,*) 'RANGULAR_MPI'
        WRITE(*,*) 'This program performs angular integration'
        WRITE(*,*) 'Input file:  rcsf.inp'
        WRITE(*,*) 'Outputfiles: mcp.30, mcp.31, .... '
        WRITE(*,*)
      END IF



*=======================================================================
*  Start mpi --- get processor info: myid, nprocs, host name; and print
*=======================================================================

      CALL startmpi2 (myid, nprocs, host, lenhost, ncount1,
     &                     startdir, permdir, tmpdir, 'MCP3_MPI')
      WRITE (idstring, '(I3.3)') myid
      print*, tmpdir , ' = tmpdir'

      lenperm = LEN_TRIM (permdir)
      lentmp = LEN_TRIM (tmpdir)

*=======================================================================
*  Get NDEF on node-0 and then send to all nodes   
*=======================================================================

      IF (myid .EQ. 0) THEN
         WRITE (istde,*)
         WRITE (istde,'(A)',ADVANCE='NO') 'Default settings?  (y/n) '
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
*    chkplt - compatibility of plant substitutions
*    setdbgmpi - debug output control parameters
*    setmc - machine- and installation- dependent constants
*    setsum - open the summary file
*    cslhmpi - load header of the csl file
*    setmcpmpi - open and check the  .mcp  files
*    strsum - append a summary of the inputs to the  .sum  file
*    factt - table of logarithms of factorials setup
*=======================================================================

      CALL chkplt ('GENMCPMPI')
      CALL setdbgmpi (DEBUG, permdir(1:lenperm) // '/genmcp.dbg')
      CALL setmc
      IF (NDEF .NE. 0 .AND. myid .EQ. 0) 
     &   CALL setsum (permdir(1:lenperm) // '/genmcp.sum')

      CALL cslhmpi (permdir(1:lenperm)//'/rcsf.inp'
     &              , ncore, nblk0, idblk)

      RESTRT = .FALSE.

      CALL setmcpmpi (myid, nprocs, ncore, idblk, 'mcp' // idstring)
      IF (NDEF .NE. 0 .AND. myid .EQ. 0) CALL strsum
      CALL factt

*=======================================================================
*     For each block, generate and sort the data
*=======================================================================

      DO nb = 1, nblock
         ncf = ncfblk(nb)        ! This ncf goes to common
         IF (myid .EQ. 0) THEN
            PRINT *
            PRINT *, 'Block ', nb, ',  ncf = ', ncf
         ENDIF

 	      !*** Load current CSL block. Memories de-allocated in mcp ***
         CALL alloc (pntriq, nnnwp  *ncf, 4)
         CALL alloc (pntjqs, nnnwp*3*ncf, 4)
         CALL alloc (pnjcup, nnnwp  *ncf, 4)

         CALL lodcslmpi (21, ncore, nb)

         !*** Open tmp.xx files for block nb ***
         CALL SETTMP (nb, kmax, 'tmp' // idstring)

         !*** Generation of MCP coefficients ***
         CALL mcpmpi (nb, RESTRT, myid, nprocs, 'mcp' // idstring)
      ENDDO

      CLOSE (24)              ! Summary file
      IF (DEBUG) CLOSE (99)   ! Debug file

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================

      CALL stopmpi2 (myid, nprocs, host, lenhost, 
     &                     ncount1, 'MCP3_MPI')

      STOP
      END
