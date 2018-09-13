************************************************************************
*                                                                      *
      SUBROUTINE PARSJL (NPEEL,ITJ,NTJ)
*                                                                      *
*   READs and parses a string that specifies J values.                 *
*                                                                      *
*   Call(s) to: [LIB92] CONVRT.                                        *
*                                                                      *
*   Written by Farid A Parpia               Last revised 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      LOGICAL NPEVEN
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 CTEGER
      CHARACTER*1 RECI
*
      include 'parameters.def'
CGG      INTEGER NNNTJV
CGG      PARAMETER (NNNTJV = 10)
      DIMENSION ITJ(NNNTJV)
*
      NPEVEN = MOD (NPEEL,2) .EQ. 0
*
      GOTO 2
*
*   Error return
*
    1 PRINT *, ' redo ...'
*
*   Read a record
*
    2 READ (*,'(A)') RECORD
*
*   Initialise NTJ
*
      NTJ = 0
*
*   Parse RECORD from left to right
*
      ISTART = 0
      I = 1
    3 RECI = RECORD(I:I)
      IF ((RECI .NE. ' ') .AND. (RECI .NE. ',')) THEN
         IF (ISTART .EQ. 0) THEN
            ISTART = I
            IFRAC = 0
         ELSE
            IF (RECI .EQ. '/') THEN
               IFRAC = I
            ENDIF
         ENDIF
      ELSE
         IF (ISTART .NE. 0) THEN
            NTJ = NTJ+1
            IF (NTJ .GT. NNNTJV) THEN
               PRINT *, 'PARSJL: Number of atomic angular'
     :               //' momenta exceeds dimensioned allowance:'
               PRINT *,' plant NTJV was set to NNNTJV.'
               STOP
            ENDIF
            IEND = I-1
            IF (IFRAC .EQ. 0) THEN
               CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS) J
               IF ((IOS .NE. 0) .OR. (J .LT. 0)) THEN
                  PRINT *, 'PARSJL: Angular momentum '
     :                     ,RECORD(ISTART:IEND),
     :                    ' incorrect;'
                  GOTO 1
               ENDIF
               J = 2*J
            ELSE
               CALL CONVRT (IEND-IFRAC,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(IFRAC+1:IEND),FMT = FORM,IOSTAT = IOS) J
               IF ((IOS .NE. 0) .OR. (J .NE. 2)) THEN
                  PRINT *, 'PARSJL: Denominator '
     :                     ,RECORD(IFRAC+1:IEND),
     :                     ' incorrect;'
                  GOTO 1
               ENDIF
               CALL CONVRT (IFRAC-ISTART,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IFRAC-1),FMT = FORM,IOSTAT = IOS) J
               IF ((IOS .NE. 0) .OR. (J .LT. 0)) THEN
                  PRINT *, 'PARSJL: Angular momentum '
     :                     ,RECORD(ISTART:IFRAC-1),
     :                     ' incorrect;'
                  GOTO 1
               ENDIF
            ENDIF
            IF (NPEVEN) THEN
               IF (MOD (J,2) .NE. 0) THEN
                  PRINT *, 'PARSJL: All angular momenta should'
     :                  //' be integers;'
                  GOTO 1
               ENDIF
            ELSE
               IF (MOD (J,2) .NE. 1) THEN
                  PRINT *, 'PARSJL: All angular momenta should'
     :                  //' be half-integers;'
                  GOTO 1
               ENDIF
            ENDIF
            ITJ(NTJ) = J
            ISTART = 0
         ENDIF
      ENDIF
*
      IF (I .LT. 256) THEN
         I = I+1
         GOTO 3
      ENDIF
*
      IF (NTJ .EQ. 0) THEN
         NTJ = -1
      ENDIF
*
      RETURN
      END
