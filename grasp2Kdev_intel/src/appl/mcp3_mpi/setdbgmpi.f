************************************************************************
*                                                                      *
      SUBROUTINE setdbgmpi (DEBUG, fullname)
      IMPLICIT NONE
*                                                                      *
*   This routine calls setdbg to open the  .dbg  file and sets the 
*   arrays that control debug printout from the GENMCP program.
*                                                                      *
*   Written by Xinghong He                  Last update: 29 Jun 1998   *
*                                                                      *
************************************************************************
*
      LOGICAL DEBUG, LDBPA, LDBPG
      INTEGER NDEF
      CHARACTER *(*) fullname

      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEFAULT/NDEF

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

      IF (myid .EQ. 0) THEN
         CALL setdbg (debug, fullname)
      ENDIF

      CALL MPI_Bcast (debug,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ldbpa,5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (ldbpg,5,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)

      RETURN
      END
