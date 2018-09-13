************************************************************************
*                                                                      *
      SUBROUTINE SETRWF
*                                                                      *
*   Open, check, load data from and close the  .rwf  file.             *
*                                                                      *
*   Call(s) to: [LIB92]: LODRWF, OPENFL.
*
*   This routine is to be replaced by setrwfa or setrwfmpi
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL FOUND
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*6 G92RWF
      CHARACTER*3 STATUS
*
*   File  .rwf  file is UNFORMATTED; it must exist
*
      DEFNAM = 'rwfn.inp'
      FORM = 'UNFORMATTED'
      STATUS = 'OLD'
*
      CALL OPENFL (23,DEFNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) THEN
         PRINT *, 'setrwf: Error opening rwfn.inp'
         STOP
      ENDIF
*
*   Check the file; if not as expected, try again
*
      READ (23,IOSTAT = IOS) G92RWF
      IF ((IOS .NE. 0) .OR. (G92RWF .NE. 'G92RWF')) THEN
         PRINT *, 'This is not a Radial WaveFunction File;'
         CLOSE (23)
         STOP
      ENDIF
*
*   Attempt to load the radial wavefunctions; if this fails,
*   try again
*
      CALL LODRWF (IERR)
      IF (IERR .NE. 0) THEN
         PRINT *, 'Radial wavefunctions defined in Configuration'
         PRINT *, ' Symmetry List File not found in Radial'
         PRINT *, ' WaveFunction File;'
         CLOSE (23)
         STOP
      ENDIF
*
*   Close the  .rwf  file
*
      CLOSE (23)
*
      RETURN
      END
