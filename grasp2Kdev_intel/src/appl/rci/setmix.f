************************************************************************
*                                                                      *
      SUBROUTINE setmix (name, idblk)
*                                                                      *
*   Opens the  .mix  file on stream 25; writes a header to this file;  *
*   calls LODMIX to interactively determine the eigenpairs required.   *
*                                                                      *
*   Call(s) to: [LIB92]: OPENFL.                                       *
*               [RCI92]: LODMIX.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 18 Dec 1992   *
*   Modified by Xinghong He               Last revision: 23 Jun 1998   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      CHARACTER*(*) name, idblk(*)*8
      CHARACTER(LEN=*), PARAMETER:: FORM = 'UNFORMATTED',
     &										STATUS = 'UNKNOWN'

!-----------------------------------------------------------------------
      k = INDEX (name,' ')
      CALL openfl (25, name(1:k-1)//'.cm', FORM, STATUS, ierr)
      IF (ierr .NE. 0) THEN
         STOP 'setmix: Error when opening .cm file'
      ENDIF

      WRITE (25) 'G92MIX'

      CALL LODMIX (idblk)

      RETURN
      END
