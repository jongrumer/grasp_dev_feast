      SUBROUTINE cofpotmpi (EOL, J, npts)
      !IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      IMPLICIT NONE

      LOGICAL EOL
      INTEGER J, npts
!Rasa start
      DOUBLE PRECISION aaa
!Rasa end

      include 'parameters.def'
CGG      INTEGER NNNP,NNN1
CGG      PARAMETER (NNNP = 590)
CGG      PARAMETER (NNN1 = NNNP+10)
      DOUBLE PRECISION YP, XP, XQ
      COMMON/POTE/YP(NNNP),XP(NNNP),XQ(NNNP)

      !POINTER (pnttmpmpi, tmpmpi(1))  ! tmp buffer for mpi operation
*P    DOUBLE PRECISION tmpmpi(npts)    ! tmp buffer for mpi operation

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      CALL SETCOF (EOL, J)
      CALL YPOT (J)
      CALL XPOT (J)
      CALL LAGCON (J, nprocs)
      CALL DACON

!Rasa start
!     print*,"cofpotmpi: myid=",myid," before red.YP(1)=",YP(1)
!     call flush(6)
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr); 
!     if (myid .EQ. 0) read *, aaa 
!Rasa end

      CALL MPI_Allreduce (MPI_IN_PLACE, YP, npts, MPI_DOUBLE_PRECISION,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
*P    CALL MPI_Allreduce (YP, tmpmpi, npts, MPI_DOUBLE_PRECISION,
*P   &  MPI_SUM, MPI_COMM_WORLD, ierr)
*P    CALL dcopy (npts, tmpmpi, 1, YP, 1)

!Rasa start
!     print*,"cofpotmpi: myid=",myid," after red.YP(1)=",YP(1)
!     call flush(6)
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr); 
!     if (myid .EQ. 0) read *, aaa 
!Rasa end

      CALL MPI_Allreduce (MPI_IN_PLACE, XP, npts, MPI_DOUBLE_PRECISION,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
*P    CALL MPI_Allreduce (XP, tmpmpi, npts, MPI_DOUBLE_PRECISION,
*P   &  MPI_SUM, MPI_COMM_WORLD, ierr)
*P    CALL dcopy (npts, tmpmpi, 1, XP, 1)

      CALL MPI_Allreduce (MPI_IN_PLACE, XQ, npts, MPI_DOUBLE_PRECISION,
     &  MPI_SUM, MPI_COMM_WORLD, ierr)
*P    CALL MPI_Allreduce (XQ, tmpmpi, npts, MPI_DOUBLE_PRECISION,
*P   &  MPI_SUM, MPI_COMM_WORLD, ierr)
*P    CALL dcopy (npts, tmpmpi, 1, XQ, 1)

      RETURN
      END



