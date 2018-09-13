************************************************************************
*                                                                      *
      SUBROUTINE PRSRSL (NORB)
*                                                                      *
*   READs and parses a list of relativistic subshell labels delimit-   *
*   ed either by blanks or by commas. All subhells associated with a   *
*   given principal quantum number can be indicated by the  wildcard   *
*   character  * .                                                     *
*                                                                      *
*   Call(s) to: DECRSL.                                                *
*               [LIB92] CONVRT.                                        *
*                                                                      *
*   Written by Farid A. Parpia             Last revised: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 CTEGER,SYM,SYMLST
      CHARACTER*1 RECI
*
      DIMENSION SYMLST(16)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/ORBNUM/NP(NNNW),N2J(NNNW),NL(NNNW)
     :      /ORBSYM/SYM(NNNW)
*
      DATA SYMLST/'s ','p-','p ','d-','d ','f-','f ','g-','g ',
     :                 'h-','h ','i-','i ','k-','k ','l-'/
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
            IF (RECI .EQ. '*') THEN
               IEND = IEND-1
               CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS) NPRIN
               IF (IOS .NE. 0) THEN
                  PRINT *, 'PRSRSL: Principal quantum number '
     :                     ,RECORD(ISTART:IEND),' could not'
     :                  //' be decoded;'
                  GOTO 1
               ENDIF
               LOC = NORB
               NORB = NORB+2*NPRIN-1
               IF (NORB .GT. NNNW) THEN
                  PRINT *, 'PRSRSL: Number of subshells exceeds'
     :                  //' dimensioned allowance: plant NW was'
     :                  //' set to NNNW.'
                  STOP
               ENDIF
               DO 4 II = 1,2*NPRIN-1
                  LOC = LOC+1
                  NP(LOC) = NPRIN
                  SYM(LOC) = SYMLST(II)
                  CALL DECRSL(LOC,IERR)
                  IF (IERR .NE. 0) THEN
                     PRINT *, 'PRSRSL: Weird error.'
                     STOP
                  ENDIF
    4          CONTINUE
            ELSE
               NORB = NORB+1
               IF (NORB .GT. NNNW) THEN
                  PRINT *, 'PRSRSL: Number of subshells exceeds'
     :                  //' dimensioned allowance: plant NW was'
     :                  //' set to NNNW.'
                  STOP
               ENDIF
               IF (RECI .EQ. '-') THEN
                  READ (RECORD(IEND-1:IEND),'(1A2)') SYM(NORB)
                  IEND = IEND-2
               ELSE
                  SYM(NORB)(2:2) = ' '
                  READ (RECI,'(1A1)') SYM(NORB)(1:1)
                  IEND = IEND-1
               ENDIF
               CALL DECRSL (NORB,IERR)
               IF (IERR .EQ. 1) THEN
                  GOTO 1
               ENDIF
               CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS)
     :            NP(NORB)
               IF (IOS .NE. 0) THEN
                  PRINT *, 'PRSRSL: Principal quantum number '
     :                     ,RECORD(ISTART:IEND),' could not'
     :                  //' be decoded;'
                  GOTO 1
               ENDIF
               IF (NL(NORB)+1 .GT. NP(NORB)) THEN
                  PRINT *, 'PRSRSL: Incorrect subshell label;'
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
