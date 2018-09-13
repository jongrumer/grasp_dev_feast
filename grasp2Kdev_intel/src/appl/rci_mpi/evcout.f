************************************************************************
*                                                                      *
      SUBROUTINE EVCOUT
*                                                                      *
*   Routine for printing eigenvectors.                                 *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*                                                                      *
*   Written by Farid A Parpia             Last revision: 06 Oct 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8          (A-H, O-Z)
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
*
      POINTER (PNILOC,ILOC(*))
*
      POINTER (PNEVEC,EVEC(*))
      POINTER (PNIVEC,IVEC(*))
*
      COMMON/EIGVEC/PNEVEC
     :      /ORB2/NCF,NW,PNTRIQ
     :      /PRNT/NVEC,PNIVEC,NVECMX
*
*   Header
*
      WRITE (24,300)
*
*   Determine which eigenvectors are to be printed out
*
      NVEX = NVEC
      CALL ALLOC (PNILOC,NVEC,4)
      DO 1 I = 1,NVEC
         ILOC(I) = IVEC(I)
    1 CONTINUE
*
*   There are eight columns across the width of a page
*
      JP = 8
*
*   Set up for the first set of columns
*
      J1 = 1
      J2 = JP
*
*   Loop over sets of columns
*
    3 IF (J2 .GT. NVEX) J2 = NVEX
*
*   Write out the column numbers; skip a line
*
      WRITE (24,301) (ILOC(L),L = J1,J2)
      WRITE (24,302)
*
*   Write out the rows
*
      DO 4 M = 1,NCF
*
         WRITE (24,303) (EVEC(M+(L-1)*NCF),L = J1,J2)
*
*   Skip a line after every tenth row
*
         IF (MOD(M,10) .EQ. 0) WRITE (24,302)
*
    4 CONTINUE
*
*   Set up for the next set of columns
*
      WRITE (24,303)
      IF (J2 .LT. NVEX) THEN
         J1 = J1+JP
         J2 = J2+JP
         GOTO 3
      ENDIF
*
*   Deallocate storage for ILOC
*
      CALL DALLOC (PNILOC)
*
      RETURN
*
  300 FORMAT (//' Eigenvectors:'/)
  301 FORMAT (1X,8(I8,7X))
  302 FORMAT (1X)
  303 FORMAT (1X,1P,8D15.7)
*
      END
