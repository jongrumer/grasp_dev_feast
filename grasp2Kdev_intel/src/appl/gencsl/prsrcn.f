************************************************************************
*                                                                      *
      SUBROUTINE PRSRCN (NP,SYM,NCORE,NORB,IOCCS)
*                                                                      *
*   READs and parses a string that specifies a relativistic configu-   *
*   ration.                                                            *
*                                                                      *
*   Call(s) to: [LIB92] CONVRT.                                        *
*                                                                      *
*   Written by Farid A Parpia               Last revised 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 CTEGER,SYM,SYMI
      CHARACTER*1 RECI
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION NP(NNNW)
      DIMENSION SYM(NNNW)
      DIMENSION IOCCS(NNNW)
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
*   Initialise IOCCS
*
      DO 3 I = 1,NORB
         IOCCS(I) = 0
    3 CONTINUE
*
*   Parse RECORD from left to right
*
      ISTART = 0
      I = 1
    4 RECI = RECORD(I:I)
      IF     (RECI .EQ. '(') THEN
         IEND = I-1
         RECI = RECORD(IEND:IEND)
         IF (RECI .EQ. '-') THEN
            READ (RECORD(IEND-1:IEND),'(1A2)') SYMI
            IEND = IEND-2
         ELSE
            SYMI(2:2) = ' '
            READ (RECI,'(1A1)') SYMI(1:1)
            IEND = IEND-1
         ENDIF
         CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
         FORM = '(1I'//CTEGER(1:LTEGER)//')'
         READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS) NPI
         IF (IOS .NE. 0) THEN
            PRINT *, 'PRSRCN: Principal quantum number '
     :               ,RECORD(ISTART:IEND),
     :              ' could not be decoded;'
            GOTO 1
         ENDIF
         DO 5 J = NCORE+1,NORB
            IF ((NP(J) .EQ. NPI) .AND. (SYM(J) .EQ. SYMI)) THEN
               ISHELL = J
               GOTO 6
            ENDIF
    5    CONTINUE
         PRINT *, 'PRSRCN: Subshell not in peel list;'
         GOTO 1
    6    IOSTRT = I+1
      ELSEIF (RECI .EQ. ')') THEN
         IOEND = I-1
         IF (IOEND .GE. IOSTRT) THEN
            CALL CONVRT (IOEND-IOSTRT+1,CTEGER,LTEGER)
            FORM = '(1I'//CTEGER(1:LTEGER)//')'
            READ (RECORD(IOSTRT:IOEND),FMT = FORM) IOCCS(ISHELL)
         ELSE
            IOCCS(ISHELL) = 0
         ENDIF
         ISTART = 0
      ELSE
         IF (ISTART .EQ. 0) ISTART = I
      ENDIF
*
      IF (I .LT. 256) THEN
         I = I+1
         GOTO 4
      ENDIF
*
      RETURN
      END
