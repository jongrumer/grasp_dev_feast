************************************************************************
      SUBROUTINE SETRWFmpi (NAME)
      !IMPLICIT NONE   ! T3E has problem reading IBM files
      IMPLICIT REAL*8                  (A-H, O-Z)
      CHARACTER NAME*(*)

*   Open, check, load data from and close the  .rwf  file.             *
*
*   Used by rcimpivu, rscfmpi, rcimpi
*                                                                      *
*   Call(s) to: [LIB92]: LODRWFmpi, OPENFL.                            *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*   MPI version by Xinghong He            Last revision: 06 Aug 1998   *
*
************************************************************************

      INTEGER myid, nprocs, ierr, istdi, istdo, istde
      COMMON /mpi/ myid, nprocs, ierr
      COMMON/iounit/istdi,istdo,istde

      INTEGER ios, ierror
      CHARACTER G92RWF*6
!-----------------------------------------------------------------------
      IF (myid .EQ. 0) THEN 
         CALL openfl (23, name, 'UNFORMATTED', 'OLD', ierror)
         IF (ierror .EQ. 1) THEN
            WRITE (istde,*) 'Error opening', name(1:LEN_TRIM (name))
            STOP
         ENDIF

*   Check the file; if not as expected, stop.

         READ (23,IOSTAT = IOS) G92RWF
         IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
            WRITE (istde,*) 'This is not a Radial WaveFunction File;'
            CLOSE (23)
            STOP
         ENDIF
      ENDIF

*   Attempt to load the radial wavefunctions

      CALL LODRWFmpi (ierror)

      IF (ierror .NE. 0) THEN
         IF (myid .EQ. 0) THEN
            WRITE (istde,*) 'Radial wavefunctions defined in CSL file'
     &                    , ' not found in Radial WaveFunction File'
            CLOSE (23)
         ENDIF
         STOP
      ENDIF

      IF (myid .EQ. 0) CLOSE (23)

      RETURN
      END
