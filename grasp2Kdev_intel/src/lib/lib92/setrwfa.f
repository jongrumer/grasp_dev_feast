************************************************************************
*                                                                      *
      SUBROUTINE setrwfa (name)
      !IMPLICIT NONE   !T3E has problem reading IBM files
      IMPLICIT REAL*8                    (A-H, O-Z)
      CHARACTER name*(*)
*                                                                      *
*   Open, check, load data from and close the  .rwf  file.             *
*                                                                      *
*   Call(s) to: [LIB92]: LODRWF, OPENFL.                               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*   Modified    by Xinghong He            Last revision: 09 Jul 1998   *
*                                                                      *
************************************************************************
*
      CHARACTER G92RWF*6
      INTEGER ierr, ios
      INTEGER istdi, istdo, istde
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------
C      print *,name,'setrwfa'
      CALL openfl (23, name, 'UNFORMATTED', 'OLD', ierr)
      IF (ierr .EQ. 1) THEN
         WRITE (istde,*) 'Error when opening', name(1:LEN_TRIM (name))
         STOP
      ENDIF
*
*   Check the file; if not as expected, stop.
*
      READ (23, IOSTAT = IOS) G92RWF
      IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
         WRITE (istde,*) 'This is not a Radial WaveFunction File;'
         CLOSE (23)
         STOP
      ENDIF
*
*   Attempt to load the radial wavefunctions; if this fails, stop
*
      CALL LODRWF (IERR)

      IF (IERR .NE. 0) THEN
         WRITE (istde,*) 'Radial wavefunctions defined in CSL file'
     &                 , ' not found in Radial WaveFunction File'
         CLOSE (23)
         STOP
      ENDIF

      CLOSE (23)

      RETURN
      END
