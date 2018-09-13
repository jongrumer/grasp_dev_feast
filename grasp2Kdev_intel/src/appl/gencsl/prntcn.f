************************************************************************
*                                                                      *
      SUBROUTINE PRNTCN (NP,SYM,IOCCS,NCORE,NORB)
*                                                                      *
*   Prints configurations to default unit.                             *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 SYM
      CHARACTER*1 CNUM
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION NP(NNNW),SYM(NNNW),IOCCS(NNNW),CNUM(2)
*
      DATA CNUM /'1','2'/
*
      IEND = 0
      DO 1 I = NCORE+1,NORB
         IF (IOCCS(I) .GT. 0) THEN
            IBEG = IEND+1
            LENTH = INT (LOG10 (REAL (NP(I))))+1
            IEND = IBEG+LENTH-1
            FORM = '(1I'//CNUM(LENTH)//')'
            WRITE (RECORD(IBEG:IEND),FORM) NP(I)
            IBEG = IEND+1
            IF (SYM(I)(2:2) .EQ. ' ') THEN
               IEND = IBEG
               RECORD(IBEG:IEND) = SYM(I)(1:1)
            ELSE
               IEND = IBEG+1
               RECORD(IBEG:IEND) = SYM(I)(1:2)
            ENDIF
            IBEG = IEND+1
            IEND = IBEG
            RECORD(IBEG:IEND) = '('
            IBEG = IEND+1
            LENTH = INT (LOG10 (REAL(IOCCS(I))))+1
            IEND = IBEG+LENTH-1
            FORM = '(1I'//CNUM(LENTH)//')'
            WRITE (RECORD(IBEG:IEND),FORM) IOCCS(I)
            IBEG = IEND+1
            IEND = IBEG
            RECORD(IBEG:IEND) = ')'
         ENDIF
    1 CONTINUE
*
      WRITE (*,'(A)') RECORD(1:IEND)
*
      RETURN
      END
