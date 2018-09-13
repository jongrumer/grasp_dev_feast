************************************************************************

      SUBROUTINE SETMCP (ncore, nblkin, idblk, filehead)
      IMPLICIT REAL*8          (A-H,O-Z)

*   Open, read, check the header of all .mcp files. Info for each
*   block is not accessed here.
*   File 30 stores the structure of H(DC) ; 
*   file 31 stores the  T  coefficients;  files 32, 33, ...,
*   store V(0), V(1), ...
*
*   This version works in both serial/parallel(mpi) environment. The
*   difference comes from parameters filehead (in the arg list), myid
*   and nprocs (in common/mpi/).
*
*   The following items are read from mcp files. myid and nprocs
*   read from the files are checked against the ones in COMMON/mpi/
*
*       ncore, nblock, kmaxf	            mcp.30 only
*       ncfblk(i), i=1, nblock)	         mcp.30 only
*       idblk(i), i=1, nblock)	            mcp.30 only
*
*       MCPLAB, nblock, myid, nprocs	   	all mcp files
*       nelec, ncf, nw			            all mcp files
*       DIAG, ICCUT, LFORDR	               all mcp files
*
*   Call(s) to: [LIB92]: CONVRT, OPENFL.
*
*   Written by Farid A. Parpia            Last revision: 19 Dec 1992   *
*   Modified by Xinghong He               Last revision: 06 Aug 1998   *
*
************************************************************************

      CHARACTER idblk(*)*8, filehead*(*)
      CHARACTER(LEN=LEN_TRIM (filehead)+3):: filnam

      LOGICAL   DIAG, LFORDR, FOUND, FOUND1
      CHARACTER CK*2, MCPLAB*3
*
      POINTER (PNTRIQ,RIQDUM)
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /FOPARM/ICCUT
     :      /MCPA/KMAXF
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ

      COMMON/iounit/istdi,istdo,istde
      COMMON /mpi/ myid, nprocs, ierr

      POINTER (pncfblk,ncfblk(0:*))
      COMMON/HBLOCK/nblock, pncfblk
!-----------------------------------------------------------------------
      LENFH = LEN_TRIM (filehead)

      filnam = filehead(1:LENFH) // '.30'
      OPEN (30, FILE = filnam, FORM = 'UNFORMATTED', STATUS = 'OLD',
     &                         IOSTAT = ierror)

!  Parameter ierror carries through mcp.30,... mcp.kmax

      READ (30, IOSTAT = ios) ncore, nblock, kmaxf
      ierror = ierror + ABS (ios)
      IF (nblock .GT. nblkin .OR. nblock .LT. 1) THEN
         WRITE (istde,*) 'setmcp: nblock = ',nblock
         STOP
      ENDIF

      CALL alloc (pncfblk, nblock+1, 4)
      ncfblk(0) = 0

      READ (30, IOSTAT = ios) (ncfblk(i), i=1, nblock)
      ierror = ierror + ABS (ios)
      READ (30, IOSTAT = ios) (idblk(i), i=1, nblock)
      ierror = ierror + ABS (ios)

!  Look for other mcp files

      FOUND = .TRUE.
      DO K = 31, 32 + KMAXF
         CALL CONVRT (K,CK,LCK)
         filnam = filehead(1:LENFH) // '.' // CK(1:2)
         INQUIRE (FILE = filnam, EXIST = FOUND1)
         FOUND = FOUND .AND. FOUND1
      ENDDO

      IF (.NOT. FOUND) THEN
         WRITE (istde,*) 'The mcp files do not exist'
         STOP
      ENDIF

!  Open the files; check file headers

      DO k = 30, 32 + KMAXF

         IF (k .NE. 30) THEN
            CALL CONVRT (k, CK, LCK)
            filnam = filehead(1:LENFH) // '.' // CK(1:2)
            CALL openfl (k, filnam, 'UNFORMATTED', 'OLD', ierror)
         ENDIF

         READ (k, IOSTAT = IOS) MCPLAB, nblock, myidd, nprocss

         ierror = ierror + ABS (ios)
         IF (myid .NE. myidd .OR. nprocs .NE. nprocss) THEN
            WRITE (istde,*) 'mcp files were generated under different'
     &              ,' processor configuration.'
            STOP
         ENDIF

         IF (MCPLAB .NE. 'MCP') THEN
            WRITE (istde,*) 'Not a sorted GRASP92 MCP File;'
            ierror = ierror + 1
         ENDIF

         READ (k, IOSTAT = IOS) NELEC, NCF, NW
         ierror = ierror + ABS (ios)
         READ (k, IOSTAT = IOS) DIAG, ICCUT, LFORDR
         ierror = ierror + ABS (ios)

         IF (ierror .NE. 0) THEN
            WRITE (istde,*) 'setmcp: Error accumulated , stopping...'
            DO I = 30, K
               CLOSE (I)
            ENDDO
            STOP 
         ENDIF
      ENDDO

      RETURN
      END
