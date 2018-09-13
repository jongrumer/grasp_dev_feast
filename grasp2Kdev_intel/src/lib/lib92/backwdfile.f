      SUBROUTINE backwdfile (nunit, nrec)
		IMPLICIT NONE
		INTEGER nunit, nrec, i

      DO i = 1, nrec
         BACKSPACE (nunit)
      ENDDO

      RETURN
      END
