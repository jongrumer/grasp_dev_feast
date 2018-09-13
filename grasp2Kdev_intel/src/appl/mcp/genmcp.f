
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
      PROGRAM GENMCP
      CALL genmcp1234567
      END
      SUBROUTINE genmcp1234567
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

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      INTEGER*4 IQA,JQSA,JCUPA

      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))

      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /STAT/PNTJQS,PNJCUP

      COMMON/MCPA/KMAX
      COMMON/iounit/istdi,istdo,istde

      COMMON /mpi/ myid, nprocs, ierr

! Things for timing
      INTEGER   ncount1

! cpath uses
      INTEGER       lenperm, lentmp
      CHARACTER*128 startdir, permdir, tmpdir

      write(*,*)
      write(*,*) 'RANGULAR'
      write(*,*) 'This program performs angular integration '
      write(*,*) 'Input file:  rcsf.inp'
      write(*,*) 'Outputfiles: mcp.30, mcp.31, ....'
      write(*,*)





!-----------------------------------------------------------------------
      myid = 0
      nprocs = 1
      CALL starttime (ncount1, 'MCP3')

*=======================================================================
*  Get NDEF   
*=======================================================================

      IF (myid .EQ. 0) THEN
         WRITE (istde,'(A)',ADVANCE='NO') ' Default settings?  (y/n) '
         YES = GETYN ()
         IF (YES) THEN
            NDEF = 0
         ELSE
            NDEF = 1
         ENDIF
      ENDIF

*=======================================================================
*
*  Checks and settings... Mostly done in backyard.
*
*    chkplt - compatibility of plant substitutions
*    setdbg - debug output control parameters
*    setmc - machine- and installation- dependent constants
*    setsum - open the summary file
*    cslh - load header of the csl file
*    setmcp - open and check the  .mcp  files
*    strsum - append a summary of the inputs to the  .sum  file
*    factt - table of logarithms of factorials setup
*=======================================================================

      CALL chkplt ('GENMCP')
      CALL setdbg (DEBUG, 'genmcp.dbg')
      CALL setmc
      IF (NDEF .NE. 0 .AND. myid .EQ. 0) 
     &   CALL setsum ('genmcp.sum')

      CALL cslh ('rcsf.inp', ncore, nblk0, idblk)

      RESTRT = .FALSE.

      CALL setmcp2 (myid, nprocs, ncore, idblk, 'mcp')
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

         CALL lodcsh2 (21, ncore, nb)

         !*** Open tmp.xx files for block nb ***
         CALL SETTMP (nb, kmax, 'tmp')

         !*** Generation of MCP coefficients ***
         CALL mcp (nb, RESTRT, myid, nprocs, 'mcp')
      ENDDO

      CLOSE (24)              ! Summary file
      IF (DEBUG) CLOSE (99)   ! Debug file

*=======================================================================
*  Execution finished; Statistics output
*=======================================================================

      CALL stoptime (ncount1, 'RANGULAR')

      STOP
      END
