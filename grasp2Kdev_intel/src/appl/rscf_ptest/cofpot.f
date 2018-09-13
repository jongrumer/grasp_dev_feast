      SUBROUTINE cofpot (EOL, J, npts)
      !IMPLICIT REAL*8           (A-H, O-Z)
      IMPLICIT NONE

      LOGICAL EOL
      INTEGER J, npts

      include 'parameters.def'
CGG      INTEGER NNNP
CGG      PARAMETER (NNNP = 590)
      REAL*8           YP, XP, XQ
      COMMON/POTE/YP(NNNP),XP(NNNP),XQ(NNNP)

      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr

!-----------------------------------------------------------------------
      CALL SETCOF (EOL, J)
      CALL YPOT (J)
      CALL XPOT (J)
      CALL LAGCON (J, nprocs)
      CALL DACON

      RETURN
      END
