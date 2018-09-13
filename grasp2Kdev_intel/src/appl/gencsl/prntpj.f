************************************************************************
*                                                                      *
      SUBROUTINE PRNTPJ (IPAR,NTJ,ITJ)
*                                                                      *
*   Prints parity and atomic angular momentum list to default unit.    *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*1 CNUM
      include 'parameters.def'
CGG      INTEGER NNNTJV
CGG      PARAMETER (NNNTJV = 10)
*
      DIMENSION ITJ(NNNTJV),CNUM(2)
*
      DATA CNUM /'1','2'/
*
      IF (IPAR .EQ. 1) THEN
         RECORD(1:13) = ' Even parity;'
         IBEG = 14
      ELSE
         RECORD(1:12) = ' Odd parity;'
         IBEG = 13
      ENDIF
*
      IF (NTJ .EQ. -1) THEN
         IEND = IBEG+17
         RECORD(IBEG:IEND) = ' all possible J''s;'
      ELSE
         IEND = IBEG+4
         RECORD(IBEG:IEND) = ' J = '
         DO 1 I = 1,NTJ
            ITJI = ITJ(I)
            IBEG = IEND+1
            IF (MOD (ITJI,2) .EQ. 0) THEN
               IF (ITJI .GT. 0) THEN
                  J = ITJI/2
                  LENTH = INT (LOG10 (REAL (J)))+1
                  IEND = IBEG+LENTH-1
                  FORM = '(1I'//CNUM(LENTH)//')'
                  WRITE (RECORD(IBEG:IEND),FORM) J
               ELSE
                  IEND = IBEG
                  RECORD(IBEG:IEND) = '0'
               ENDIF
            ELSE
               LENTH = INT (LOG10 (REAL (ITJI)))+1
               IEND = IBEG+LENTH-1
               FORM = '(1I'//CNUM(LENTH)//')'
               WRITE (RECORD(IBEG:IEND),FORM) ITJI
               IBEG = IEND+1
               IEND = IBEG+1
               RECORD(IBEG:IEND) = '/2'
            ENDIF
            IF (I .LT. NTJ) THEN
               IBEG = IEND+1
               IEND = IBEG+1
               RECORD(IBEG:IEND) = ', '
            ELSE
               IBEG = IEND+1
               IEND = IBEG
               RECORD(IBEG:IEND) = ';'
            ENDIF
    1    CONTINUE
      ENDIF
*
      WRITE (*,'(A)') RECORD(1:IEND)
*
      RETURN
      END
