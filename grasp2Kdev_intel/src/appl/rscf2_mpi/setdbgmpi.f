************************************************************************
*
      SUBROUTINE SETDBGmpi (dbgfile)
      IMPLICIT NONE
*
*   This subroutine calls setdbg to set the arrays that control debug
*   printout from the radial and angular modules of the GRASP92
*   This is done on node-0 (the master). the results are then
*   broadcasted to all other nodes.
*
*   Xinghong He  98-07-06
*
************************************************************************
      CHARACTER*(*) dbgfile
      LOGICAL LDBPA, LDBPG, LDBPR
      COMMON/DEBUGA/LDBPA(5)
     &      /DEBUGG/LDBPG(5)
     &      /DEBUGR/LDBPR(30)

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      IF (myid .EQ. 0) THEN
         CALL setdbg (dbgfile)
      ENDIF

      CALL MPI_Bcast (ldbpa, 5, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ldbpg, 5, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_Bcast (ldbpr,30, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)

      RETURN
      END
