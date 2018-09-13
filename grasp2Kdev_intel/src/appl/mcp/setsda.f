************************************************************************
*                                                                      *
      SUBROUTINE SETSDA (outsda, NNONZ, LPRINT, nb, myid, nprocs, fhead)
*                                                                      *
*   This routine examines lists                                        *
*                               (IC,IR,npos)                           *
*   to set up the array  IENDC  required by the Davidson eigensolver   *
*   of Stathopoulos and Fischer.                                       *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, DALLOC.                        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 10 Dec 1992   *
*   Modified for block interactions by C. F. Fischer        May 1997   *
*   Modified for multi-processors   by Xinghong He       03 Jul 1998   *
*
*  Currently shared by mcpblk, mcpmpi
*
************************************************************************
*
      EXTERNAL outsda     ! a subroutine call
      LOGICAL LPRINT
      CHARACTER fhead*(*)

      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB2/NCF,NW,PNTRIQ

      COMMON/iounit/istdi,istdo,istde
! Locals...
      POINTER (PIENDC,IENDC(0:*))
      POINTER (PNIROW,IROW(*))

      CHARACTER*3 MCPLAB
!-----------------------------------------------------------------------

      IF (myid .EQ. 0)
     & PRINT *, 'Analysing sparse matrix array definition file ...', 30

      READ (30) MCPLAB, mb
      IF (nb .NE. mb ) THEN
         WRITE (istde,*) 'setsda: nb = ', nb, '.NE. mb (=', mb,')'
         STOP
      ENDIF
*
*   Allocate storage for IENDC(0:NCF)
*
      CALL ALLOC (PIENDC, NCF+1, 4)
      CALL ALLOC (PNIROW, NNONZ, 4)
*
*   Analyse data on file 30; set up IENDC and IROW
* In multiprocessor environment, iendc of each node will have the 
* same length (ncf+1); but will have its own part filled. irow is
* local, and its length is determined by the local parameter nnonz.

      IEND = 0
      ICLAST = 0
      DO I = 1, NNONZ
         READ (30) IC, IROW(I), npos
         IF (IC .NE. ICLAST) THEN
            IENDC(ICLAST) = IEND
            ICLAST = IC
         ENDIF
         IEND = npos
      ENDDO
!xhh - changed to suits MPI environment as well
!      IENDC(NCF) = IEND
      IENDC(IC) = IEND
*
*   Sorting complete; rewrite to mcpXXX.30 file
*
      OPEN (29,FILE=fhead//'.30',STATUS='OLD',FORM='UNFORMATTED',
     :  IOSTAT=IERR, POSITION='APPEND')
      IF (IERR .NE. 0) THEN
         WRITE (istde, *) ' Error when opening the file mcp.30'
         STOP
      END IF

      WRITE (29) 'MCP', nb, ncf
      WRITE (29) NNONZ
      WRITE (29) (IENDC(I),I=myid+1,NCF,nprocs), (IROW(I),I=1,NNONZ)
      CLOSE (29)

      CALL outsda (lprint, nnonz, ncf, irow, iendc)
*
*   Deallocate storage
*
      CALL DALLOC (PIENDC)
      CALL DALLOC (PNIROW)

      RETURN
      END
