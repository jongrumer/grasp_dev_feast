************************************************************************
*                                                                      *
      SUBROUTINE PRSRCN (RECORD,NCORE,IOCCS,IERR)
*                                                                      *
*   READs and parses a string that specifies a configuration.          *
*                                                                      *
*   Written by Farid A Parpia              Last revised: 16 Oct 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 RECORD
      CHARACTER*5 FORM
      CHARACTER*2 NH,SYMI
      CHARACTER*1 CNUM,RECI
*
      DIMENSION CNUM(3)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
*
      DIMENSION IOCCS(NNNW)
*
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
      DATA CNUM/'1','2','3'/
*
*   Initialise IOCCS for the peel subshells
*
      DO 1 I = NCORE+1,NW
         IOCCS(I) = 0
    1 CONTINUE
*
*   Parse RECORD from left to right
*
      ISTART = 0
      I = 1
    2 RECI = RECORD(I:I)
      IF     (RECI .EQ. '(') THEN
         IEND = I-1
    3    RECI = RECORD(IEND:IEND)
         IF (RECI .EQ. ' ') THEN
            IEND = IEND-1
            GOTO 3
         ELSEIF (RECI .EQ. '-') THEN
            READ (RECORD(IEND-1:IEND),'(1A2)') SYMI
            IEND = IEND-2
         ELSE
            SYMI(2:2) = ' '
            READ (RECI,'(1A1)') SYMI(1:1)
            IEND = IEND-1
         ENDIF
         LENTH = IEND-ISTART+1
         FORM = '(1I'//CNUM(LENTH)//')'
         READ (RECORD(ISTART:IEND),FMT = FORM,IOSTAT = IOS) NPI
         IF (IOS .NE. 0) THEN
            PRINT *, 'PRSRCN: Principal quantum number ',
     :                RECORD(ISTART:IEND)
            PRINT *, ' could not be decoded.'
            IERR = 1
            GOTO 6
         ENDIF
         DO 4 J = NCORE+1,NW
            IF ((NP(J) .EQ. NPI) .AND. (NH(J) .EQ. SYMI)) THEN
               ISHELL = J
               GOTO 5
            ENDIF
    4    CONTINUE
         PRINT *, 'PRSRCL: Not a peel subshell.'
         IERR = 2
         GOTO 6
    5    IOSTRT = I+1
      ELSEIF (RECI .EQ. ')') THEN
         IOEND = I-1
         LENTH = IOEND-IOSTRT+1
         FORM = '(1I'//CNUM(LENTH)//')'
         READ (RECORD(IOSTRT:IOEND),FMT = FORM,IOSTAT = IOS) IOCCI
         IF (IOS .NE. 0) THEN
            PRINT *, 'PRSRCN: Occupation number ',
     :                RECORD(IOSTRT:IOEND)
            PRINT *, ' could not be decoded.'
            IERR = 3
            GOTO 6
         ENDIF
         NKJI = NKJ(ISHELL)
         IF (NKJI .LE. 7) THEN
            IF ((IOCCI .LT. 0) .OR. (IOCCI .GT. NKJ(ISHELL)+1)) THEN
               PRINT *, 'PRSRCN: Occupation specified'
               PRINT *, ' incorrectly for ',NP(ISHELL),NH(ISHELL)
               PRINT *, ' subshell.'
               IERR = 4
               GOTO 6
            ENDIF
         ELSE
            IF ((IOCCI .LT. 0) .OR. (IOCCI .GT. 2)) THEN
               PRINT *, 'PRSRCN: Occupation specified'
               PRINT *, ' incorrectly for ',NP(ISHELL),NH(ISHELL)
               PRINT *, ' subshell.'
               IERR = 5
               GOTO 6
            ENDIF
         ENDIF
         IOCCS(ISHELL) = IOCCI
         ISTART = 0
      ELSE
         IF (ISTART .EQ. 0) ISTART = I
      ENDIF
*
      IF (I .LT. 256) THEN
         I = I+1
         GOTO 2
      ENDIF
*
      IERR = 0
*
    6 RETURN
      END
