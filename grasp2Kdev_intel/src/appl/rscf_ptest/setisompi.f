************************************************************************
      SUBROUTINE setisompi (isofile)
      IMPLICIT REAL*8          (A-H,O-Z)
      CHARACTER*(*) isofile

* An MPI container for setiso which opens, checks, and loads isofile
* to get isotope data. Data loaded are:
*     EMN,Z,/NPAR/,/NSMDAT/
* where /.../ means whole common block.
*
* Xinghong He 98-08-06
*
************************************************************************
      INCLUDE 'parameters.def'         ! Jon Grumer (Lund, 2013)

      COMMON/DEF1/EMN,IONCTY,NELEC,Z   ! out, ioncty and nelec not used
     :      /DEF11/FMTOAU,AUMAMU       ! in
     :      /NPAR/PARM(2),NPARM        ! out
     :      /NSMDAT/SQN,DMOMNM,QMOMB   ! out
     :      /GRID/R(NNN1),RP(NNN1),RPOR(NNN1),RNT,H,HP,N ! out - Jon Grumer, Lund 2013

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      IF (myid .EQ. 0) CALL SETISO (isofile)

      CALL MPI_Bcast (Z,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (EMN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (PARM,2,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NPARM,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (SQN,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (DMOMNM,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &  ierr)
      CALL MPI_Bcast (QMOMB,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     &  ierr)
      CALL MPI_Bcast (RNT,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (H,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (HP,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (N,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      RETURN
      END
