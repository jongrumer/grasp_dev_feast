************************************************************************
*
      SUBROUTINE setcslmpi (name, ncore, idblk)
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      INTEGER   ncore
      CHARACTER idblk(*)*8, name*(*)
*
*  A container which calls lib92/lodcsh to get
*     ncore, nelec, nw, np(), nak(), nkl(), nkj(), nh()
*  It then calls lodcsl to load the entire file (all blocks), obtaining
*     IQA, JQSA, JCUPA
*  Note that unlike setcsl of mcp[mpi] and rci[mpi]vu, the following
*  quantities are assumed to be known:
*     nblock, ncfblk(), idblk(), ncftot.
*  In this rscfmpi, they have been read into memory from mcp files via 
*  a call to setmcp.
*
*  Xinghong He 98-08-06
*
************************************************************************
*
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWP = 30)
      CHARACTER str*15

      POINTER (pncfblk, ncfblk(0:*))
      COMMON/hblock/nblock, pncfblk

      INTEGER*4 IQA,JQSA,JCUPA
      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))

      CHARACTER*2 NH
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCFtot,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

! node-0 opens, reads the header of the file

      IF (myid .EQ. 0) THEN
         CALL OPENFL (21, name, 'FORMATTED', 'OLD', ierr)
         IF (ierr .EQ. 1) THEN
            PRINT *, 'Error when opening ',name(1:LEN_TRIM (name))
            STOP
         ENDIF

         READ (21,'(1A15)',IOSTAT = ios) str
         IF ((ios .NE. 0) .OR.
     :         (str .NE. 'Core subshells:')) THEN
            PRINT *, 'Not a Configuration Symmetry List File;'
            CLOSE (21)
            STOP
         ENDIF

         !..Load header of <name> file
         CALL LODCSH (21, NCORE)
      ENDIF

! Broadcast results to other nodes. ncfblk should be allocated
! on these nodes (with myid .ne. 0)

      CALL MPI_Bcast (nw,    1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ncore, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nelec, 1, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (np(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nak(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkl(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nkj(1),nw, MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (nh(1),2*nw, MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)

! Allocate memories for all blocks and then load the entire file

      CALL ALLOC (PNTRIQ,NNNWP  *ncftot,4)
      CALL ALLOC (PNTJQS,NNNWP*3*ncftot,4)
      CALL ALLOC (PNJCUP,NNNWP  *ncftot,4)

      CALL lodcslmpi (21, ncore, -119)
      ! -119 means "load all blocks"

      RETURN
      END
