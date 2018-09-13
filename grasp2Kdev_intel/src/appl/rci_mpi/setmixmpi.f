************************************************************************
*                                                                      *
      SUBROUTINE setmixmpi (name, idblk)
*                                                                      *
*   Opens the  .mix  file on stream 25; writes a header to this file;  *
*   calls LODMIXmpi to interactively determine the eigenpairs required. 
*                                                                      *
*   Call(s) to: [LIB92]: OPENFL.                                       *
*               [RCI92]: LODMIXmpi                                     *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*   MPI version by Xinghong He            Last revision:  9 May 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*(*) name, idblk(*)*8
      CHARACTER(LEN=*), PARAMETER:: FORM = 'UNFORMATTED',
     &                              STATUS = 'UNKNOWN'

      INCLUDE 'mpif.h'
      COMMON /mpi/ myid, nprocs, ierr
!-----------------------------------------------------------------------
      IF (myid .EQ. 0) THEN
         k = INDEX (name,' ')
         CALL openfl (25, name(1:k-1)//'.cm', FORM, STATUS, ierr)
         IF (ierr .NE. 0) THEN
            CALL stopmpi ('setmix: Error when opening .cm file', myid)
         ENDIF

         WRITE (25) 'G92MIX'
      ENDIF

      CALL LODMIXmpi (idblk)

      RETURN
      END
