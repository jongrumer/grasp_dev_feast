************************************************************************
      SUBROUTINE lodcslmpi (nfile, ncore, jblock)
      IMPLICIT NONE
      INTEGER nfile, ncore, jblock

* An MPI container of lodcsh2 which loads CSL list of the current block 
* into memory. It forwards the call together with the same set of 
* parameters to lodcsh2 and then broadcasts the results to all nodes.
*
* Note: Memories have been allocated/deallocated each block outside.
* This subroutine calls lodcsh2 on node-0 to generate the data for the 
* block; and then broadcasts to all other nodes. A new MPI data type
* of 4 byte-long is created to handle 64-bit machines whose MPI
* implementation does not support 4-byte integers. If jblock=-119,
* then ALL blocks will be loaded instead of just one. This is 
* implemented in lodcsh2.
*
* Currently used by rcimpivu, mcpmpi, rscfmpivu 
*
* Xinghong He 98-08-06
*
************************************************************************

CGG      INTEGER   NNNWP
      include 'parameters.def'
CGG      PARAMETER (NNNWP = 30)

      INTEGER   NCF, NW
      INTEGER*4 IQA,JQSA,JCUPA, MPIX_INT4

      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))
      COMMON/ORB2/NCF,NW,PNTRIQ
     &      /STAT/PNTJQS,PNJCUP

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

      IF (myid .EQ. 0) THEN
         CALL lodcsh2 ((nfile), (ncore), (jblock))
      ENDIF

! Construct mpi data type for Integer*4 and then broadcast.

      CALL mpix_bytes (4, MPIX_INT4, ierr)

      CALL MPI_Bcast (IQA,   NNNWP*NCF,MPIX_INT4,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (JQSA,3*NNNWP*NCF,MPIX_INT4,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (JCUPA, NNNWP*NCF,MPIX_INT4,0,MPI_COMM_WORLD,ierr)

      RETURN
      END
