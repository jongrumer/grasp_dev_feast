************************************************************************
*                                                                      *
      SUBROUTINE FNAME(NAME)
*                                                                      *
*   Determines the name of the initial and final states                *
*                                                                      *
*   Written by Per Jonsson                                             *       
*                                                                      *
************************************************************************
*
      CHARACTER*24 NAME(2)
*
*   Obtain the names of the initial and final state files
*
    1 PRINT *, ' Name of the Initial state'
      READ(*,'(A)') NAME(1)
      PRINT *, ' Name of the Final state'
      READ(*,'(A)') NAME(2)
*
      J = INDEX(NAME(1),' ')
      IF (J .EQ. 1) THEN
         PRINT *, ' Names may not start with blanks'
         GOTO 1
      ENDIF
*
      J = INDEX(NAME(2),' ')
      IF (J .EQ. 1) THEN
         PRINT *, ' Names may not start with blanks'
         GOTO 1
      ENDIF

      RETURN 
      END
