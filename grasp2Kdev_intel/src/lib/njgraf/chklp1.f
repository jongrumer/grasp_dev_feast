************************************************************************
*                                                                      *
      SUBROUTINE CHKLP1 (FAIL)
*                                                                      *
*   This routine checks if  there are active triads with two identi-   *
*   cal  arguments.  This is a loop of  order 1 (a "lollypop").  The   *
*   other argument  must then be zero;  i.e. j1(j) = 1 in 2j+1 nota-   *
*   tion. Suppression of the loop introduces factors and phases. Two   *
*   Two triads  become inactive.  All this is  performed by invoking   *
*   ZERO with first argument 1.                                        *
*                                                                      *
*   Written by  Marcel Klapisch for correcting  an error detected by   *
*   Charlotte F Fischer.    This version includes Marcel's fix of 12   *
*   June 1992.                                                         *
*                                         Last revision: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
      INTEGER ARROW
      LOGICAL FAIL,FREE,TABS
      CHARACTER*6 NAME,NAMSUB
*
      PARAMETER (MANGM = 60,MTRIAD = 12,M2TRD = 2*MTRIAD)
*
      COMMON/COUPLE/M,N,J1(MANGM),J2(MTRIAD,3),J3(MTRIAD,3),
     :              FREE(MANGM)
     :      /DEBUG/IBUG1,IBUG2,IBUG3,IBUG4,IBUG5,IBUG6
     :      /NAM/NAMSUB
     :      /TREE/J23(M2TRD,3),ARROW(M2TRD,3),LINE(MANGM,2),
     :            LCOL(MANGM,2),TABS(M2TRD),NBTR
*
      DATA NAME/'CHKLP1'/
*
      NAMSUB = NAME
      NBTR1 = 2*(N-1)
      DO 10 L = 1,NBTR1
      IF (.NOT. TABS(L)) THEN
         JDIF = 0
         IF     (J23(L,1) .EQ. J23(L,2)) THEN
            JDIF = J23(L,3)
         ELSEIF (J23(L,1) .EQ. J23(L,3)) THEN
            JDIF = J23(L,2)
         ELSEIF (J23(L,2) .EQ. J23(L,3)) THEN
            JDIF = J23(L,1)
         ENDIF
         IF (JDIF .NE. 0) THEN
*
*   Putting the link to 0. ZERO changes NBTR
*
            FAIL = .FALSE.
            IF (J1(JDIF) .NE. 1.AND. .NOT. FREE(JDIF)) THEN
               FAIL = .TRUE.
               IF (IBUG3 .EQ. 1) WRITE (99,300) JDIF,J1(JDIF)
               RETURN
            ELSE
               CALL ZERO (1,JDIF,FAIL)
               IF (FAIL) RETURN
            ENDIF
         ENDIF
      ENDIF
  10  CONTINUE
      IF (JDIF .NE. 0) CALL PRINTJ (NAME,4)
*
      RETURN
*
  300 FORMAT (1X,'JDIF = ',1I2,'; should be 0; J1(JDIF) = ',1I2,
     :           '; RECUP -> 0.')
*
      END
