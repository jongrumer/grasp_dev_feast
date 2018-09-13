************************************************************************
*                                                                      *
      SUBROUTINE NEIBOR (LC,L1,L2)
*                                                                      *
*   Gives the positions of the other two arguments in the triad.       *
*                                                                      *
*                                           Last update: 16 Oct 1992   *
*                                                                      *
************************************************************************
      IMPLICIT REAL*8          (A-H, O-Z)
*
      IF (LC .LT. 2) THEN
         L1 = 2
         L2 = 3
      ELSEIF (LC .EQ. 2) THEN
         L1 = 3
         L2 = 1
      ELSE
         L1 = 1
         L2 = 2
      ENDIF
*
      RETURN
      END
