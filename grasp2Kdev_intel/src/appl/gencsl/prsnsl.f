************************************************************************
*                                                                      *
      SUBROUTINE PRSNSL (NORB)
*                                                                      *
*   READs and parses a list of  nonrelativistic subshell  labels de-   *
*   limited either by blanks or commas. The wildcard  *  can be used   *
*   to abbreviate all subshells for a given principal quantum number.  *
*                                                                      *
*   Call(s) to: CONVRT.                                                *
*                                                                      *
*   Written by Farid A. Parpia             Last revised: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 CTEGER
      CHARACTER*1 NRSYM,RECI,SYMLST
*
      DIMENSION SYMLST(9)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/NRORBN/NPNR(NNNW),N2LNR(NNNW)
     :      /NRORBS/NRSYM(NNNW)
*
      DATA SYMLST/'s','p','d','f','g','h','i','k','l'/
*
      NORSAV = NORB
      GOTO 2
*
*   Error return
*
    1 PRINT *, ' redo ...'
      NORB = NORSAV
*
*   Read a record
*
    2 READ (*,'(A)') RECORD
*
*   Parse RECORD from left to right
*
      ISTART = 0
      I = 1
    3 RECI = RECORD(I:I)
      IF ((RECI .NE. ' ') .AND. (RECI .NE. ',')) THEN
         IF (ISTART .EQ. 0) ISTART = I
      ELSE
         IF (ISTART .NE. 0) THEN
            IEND = I-1
            RECI = RECORD(IEND:IEND)
            IEND = IEND-1
            CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
            FORM = '(1I'//CTEGER(1:LTEGER)//')'
            READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS) NP
            IF (IOS .NE. 0) THEN
               PRINT *, 'PRSNSL: Principal quantum number '
     :                  ,RECORD(ISTART:IEND),' could not'
     :               //' be decoded;'
               GOTO 1
            ENDIF
            IF (RECI .EQ. '*') THEN
               NORB = NORB+NP
               IF (NORB .GT. NNNW) THEN
                  PRINT *, 'PRSNSL: Number of subshells exceeds'
     :                  //' dimensioned allowance: plant NW was'
     :                  //' set to NNNW.'
                  STOP
               ENDIF
               LOC = NORB-NP
               DO 4 II = 1,NP
                  LOC = LOC+1
                  NPNR(LOC) = NP
                  NRSYM(LOC) = SYMLST(II)
                  CALL DECNSL (LOC,IERR)
                  IF (IERR .EQ. 1) THEN
                     PRINT *, 'PRSNSL: Weird error.'
                     STOP
                  ENDIF
    4          CONTINUE
            ELSE
               NORB = NORB+1
               IF (NORB .GT. NNNW) THEN
                  PRINT *, 'PRSNSL: Number of subshells exceeds'
     :                  //' dimensioned allowance: plant NW'
     :                  //' was set to NNNW.'
                  STOP
               ENDIF
               NPNR(NORB) = NP
               NRSYM(NORB) = RECI
               CALL DECNSL (NORB,IERR)
               IF (IERR .EQ. 1) THEN
                  GOTO 1
               ENDIF
               IF (N2LNR(NORB)/2+1 .GT. NPNR(NORB)) THEN
                  PRINT *, 'PRSNSL: Incorrect subshell label;'
                  GOTO 1
               ENDIF
            ENDIF
            ISTART = 0
         ENDIF
      ENDIF
*
      IF (I .LT. 256) THEN
         I = I+1
         GOTO 3
      ENDIF
*
      RETURN
      END
