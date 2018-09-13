************************************************************************
*
      SUBROUTINE cslhmpi (name, ncore, nblkin, idblk)
      IMPLICIT REAL*8          (A-H,O-Z)
*
*  A container which calls setcsll to open, read <name> file to get 
*     nblock, ncfblk(), idblk(), ncftot. 
*  It then calls lib92/lodcsh to get
*     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
*  The file pointer points to the first CSL record after this routine.
*
*  Called by rcimpivu, mcpmpi
*
*  Xinghong He 98-06-23
*
************************************************************************
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
      CHARACTER*8 idblk(*), name*(*)

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      INTEGER*4 IQAdum
      POINTER (PNTRIQ,IQAdum)
      CHARACTER*2 NH
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCFtot,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

! node-0 does exactly the same as the serial code does

      IF (myid .EQ. 0) THEN
         CALL ALLOC   (pncfblk, nblkin+1, 4)
         CALL SETCSLL (21, name, nblkin, nblock, ncfblk(1), ncftot, 
     &                 idblk)
         CALL RALLOC  (pncfblk, nblkin+1, nblock+1, 4)
         REWIND (21)
         READ   (21,*)
         !..Load header of <name> file
         CALL LODCSH (21, NCORE)
      ENDIF

! Broadcast results to other nodes. ncfblk should be allocated
! on these nodes (with myid .ne. 0)

      CALL MPI_Bcast (nblock, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ncftot, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      IF (myid .NE. 0) THEN
      	CALL alloc (pncfblk, 1+nblock,4)
      ENDIF

      CALL MPI_Bcast (ncfblk, 1+nblock, MPI_INTEGER,0,
     &                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (idblk, 8*nblock, MPI_CHARACTER,0,
     &                                        MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ncore, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nelec, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nw,    1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (np,   nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nak,  nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkl,  nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkj,  nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nh, 2*nw, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

      RETURN
      END
