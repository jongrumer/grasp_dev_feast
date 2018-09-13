************************************************************************
      SUBROUTINE setmcpmpi (myid, nprocs, ncore, idblk, filehead)
      IMPLICIT REAL*8          (A-H,O-Z)
*
* A wrapper for setmcp/getinf. setmcp/getinf are then shared by serial
* and MPI programs.
*
*   Written by Xinghong He                Last revision: 30 Jun 1998
*
************************************************************************
      INTEGER   myid, nprocs, ncore
      CHARACTER idblk(*)*8, filehead*(*)

      POINTER (PNTRIQ,RIQDUMMY)
      LOGICAL DIAG, LFORDR
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /FOPARM/ICCUT(100)
     :      /MCPA/KMAX
     :      /MCPB/DIAG,LFORDR
     :      /ORB2/NCF,NW,PNTRIQ

      INCLUDE 'mpif.h'
      !COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

      CALL SETMCP (myid, nprocs, ncore, idblk, filehead)

! DIAG, ICCUT, LFORDR are set in GETINF

      IF (myid .EQ. 0) THEN
         CALL GETINF
      ENDIF

      CALL MPI_Bcast (ICCUT(1), 1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (DIAG,  1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (LFORDR,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      DO K = 30, 32+KMAX
         WRITE (K) NELEC,NCF,NW
         WRITE (K) DIAG,ICCUT(1),LFORDR
      ENDDO

      DO K = 30, 32+KMAX
         CLOSE (K)
      ENDDO

      RETURN
      END
