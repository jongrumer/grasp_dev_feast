************************************************************************
      SUBROUTINE posfile (mode, nunit, nrec)
      IMPLICIT NONE

      INTEGER mode, nunit, nrec, i
      INTEGER, PARAMETER:: FORWD = 0, BACKWD = 1
!----------------------------------------------------------------------

      SELECT CASE (mode)

      CASE (FORWD)
         REWIND (nunit)
         DO i = 1, nrec
            READ (nunit)
         ENDDO

      CASE (BACKWD)
         DO i = 1, nrec
            BACKSPACE (nunit)
         ENDDO

      ENDSELECT

      RETURN
      END
