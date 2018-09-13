************************************************************************
      SUBROUTINE lodcslmpijb (NAME,ncore)
c     IMPLICIT NONE
c     INTEGER NAME, ncore

* An MPI container of SETCSLAjb
* into memory. It forwards the call together with the same set of 
* parameters to SETCSLAjb and then broadcasts the results to all nodes.
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
cb
cb adapted for  hfs92MPI by Jacek Bieron
cb 14 Apr 2008
cb
*
************************************************************************

      IMPLICIT REAL*8          (A-H,O-Z)

      CHARACTER*24 NAME

c      INTEGER   NNNWP
c      PARAMETER (NNNWP = 30)
c
cbieron include 'parameters.def'
c
c      PARAMETER (NNNP = 390)
c      PARAMETER (NNN1 = 400)
c      PARAMETER (NNNW = 120)
c      
      include 'parameters.def'
c     
      INTEGER   NCF, NW
      INTEGER*4 IQA,JQSA,JCUPA, MPIX_INT4

      CHARACTER*2 NH

      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     &      /STAT/PNTJQS,PNJCUP

      INCLUDE 'mpif.h'
      INTEGER myid, nprocs, ierr
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------

      if (myid .eq. 0) then
         CALL SETCSLAjb (NAME,NCORE)
c     print *,' lodcslmpijb: myid = ',myid,' NNNWP, NCF = ',NNNWP, NCF
c     print *,' lodcslmpijb: myid = ',myid,' NW , NCF = ', NW, NCF
      endif

      call MPI_Barrier (MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast (NCF,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NW ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      CALL MPI_Bcast (NP, NW,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NAK,NW,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NKL,NW,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NKJ,NW,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (NH ,NW,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if (myid .ne. 0) then
c     print *,' lodcslmpijb: myid = ',myid,' NW , NCF = ', NW, NCF
      CALL ALLOC (PNTJQS,NNNWP*3*NCF,4)
      CALL ALLOC (PNJCUP,NNNWP  *NCF,4)
      CALL ALLOC (PNTRIQ,NNNWP  *NCF,4)
      endif

! Construct mpi data type for Integer*4 and then broadcast.

      CALL mpix_bytes (4, MPIX_INT4, ierr)

c     print *, ' after mpix_bytes '

      CALL MPI_Bcast (IQA,   NNNWP*NCF,MPIX_INT4,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (JQSA,3*NNNWP*NCF,MPIX_INT4,0,MPI_COMM_WORLD,ierr)
      CALL MPI_Bcast (JCUPA, NNNWP*NCF,MPIX_INT4,0,MPI_COMM_WORLD,ierr)

c     print *, ' after JCUPA  '
      RETURN
      END
