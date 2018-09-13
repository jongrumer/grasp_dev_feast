      SUBROUTINE outsdampi (lprint, nnonz, ncf, irow, iendc)
      IMPLICIT NONE

      LOGICAL lprint
      INTEGER nnonz, ncf, irow(nnonz), iendc(0:ncf)

      INTEGER i, ic, nnonz_a, myid, nprocs, ierr

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr

		CHARACTER*(*), PARAMETER:: myname = 'outsdampi'
!-----------------------------------------------------------------------

      CALL MPI_Reduce (nnonz, nnonz_a, 1, MPI_INTEGER, MPI_SUM, 0,
     &         MPI_COMM_WORLD, ierr)

      IF (myid .EQ. 0)
     & PRINT *, ' ... complete; density of non-zero elements of H(DC): '
     &       , nnonz_a, '/', (NCF*(NCF+1))/2
*
*   Debug printout
*
      IF (lprint) THEN

         IF (myid .EQ. 0) THEN
            WRITE (99,*)
            WRITE (99,*) 'From ', myname, ' :'
            WRITE (99,301) NNONZ_a

            PRINT *, 'This part not finished. See ', myname
            WRITE (99,*) 'This part not finished. See ', myname
         ENDIF

!         DO IC = 1, NCF
!            DO I = IENDC(IC-1) + 1, IENDC(IC)
!               WRITE (99,302) IC, IROW(I), I
!            ENDDO
!        ENDDO
      ENDIF

  301 FORMAT (' Number of nonzero elements in H(DC): ',1I4)
  302 FORMAT (' Column ',1I2,', row ',1I2,', sparse matrix index ',1I4)

      RETURN
      END
