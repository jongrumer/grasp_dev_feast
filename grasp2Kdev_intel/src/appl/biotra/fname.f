************************************************************************
*                                                                      *
      SUBROUTINE FNAME(NAME)
*                                                                      *
*   Determines the name of the initial and final states                *
*   In addition this subroutine determines which J symmetries          *
*   that are to be transformed                                         *
*                                                                      *
*   Written by Per Jonsson                                             *       
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      CHARACTER*24 NAME(2)
      LOGICAL GETYN, YES

      COMMON/JQJC/JQJ1,JQJ2,ITJQJ1(40),ITJQJ2(40),NTRANS
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

* Per april 2007
* Check if the initial and final states are identical.

      IF (TRIM(NAME(1)).EQ.TRIM(NAME(2))) THEN
         PRINT *
         PRINT *, ' Initial and final states are identical and there is'
         PRINT *, ' no need for the biorthogonal transformation. Just  '
         PRINT *, ' copy name.w to name.bw and name.(c)m to name.(c)bm '
         PRINT *, ' and run bioscl                                     '
         PRINT *


         PRINT *, ' Do you want to continue anyway ? '
         YES = GETYN ()
         IF (YES.EQV..FALSE.) STOP
      END IF

* end Per 2007


      PRINT *, ' Transformation of all J symmetries?'
      YES = GETYN () 
      IF (YES) THEN
         NTRANS = 0
      ELSE
         NTRANS = 1
         PRINT *, 
     :   ' Number of initial state J symmetries to be transformed'
         READ(*,*) JQJ1 
         PRINT *, ' Give the J symmetries in the form 2*J'
         READ(*,*) (ITJQJ1(I),I=1,JQJ1)
         DO I =1,JQJ1
            ITJQJ1(I) = ITJQJ1(I) + 1
         ENDDO
         PRINT *, 
     :   ' Number of final state J symmetries to be transformed'
         READ(*,*) JQJ2
         PRINT *, ' Give the J symmetries in the form 2*J'
         READ(*,*) (ITJQJ2(I),I=1,JQJ2)
         DO I =1,JQJ2
           ITJQJ2(I) = ITJQJ2(I) + 1
         ENDDO
      ENDIF

      RETURN 
      END
