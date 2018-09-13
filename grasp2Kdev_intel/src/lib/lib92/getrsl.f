************************************************************************
*                                                                      *
      SUBROUTINE GETRSL (indx, NSUBS)
*                                                                      *
*   READs and parses a list of relativistic subshell labels delimit-   *
*   ed either by blanks or by commas. An asterisk may be used as the   *
*   `wildcard' character.                                              *
*
*   Output:
*     NSUBS - the # of orbitals parsed
*     indx(1:NSUBS) - indeces of these orbitals
*                                                                      *
*   Call(s) to: [LIB92]: LDIGIT.                                       *
*                                                                      *
*   Written by Farid A. Parpia             Last revised: 18 Dec 1992   *
*   Modified by Xinghong He                Last revised: 09 Jul 1998   *
*                                                                      *
************************************************************************
*
      INTEGER indx(*)

      LOGICAL FOUND,LDIGIT,NLANY,NLOK,NPANY,NPOK,NSANY,NSOK
      CHARACTER*1 NLQ, NSQ, RECI, CNUM*2, NH*2, FORM*5, RECORD*500

      include 'parameters.def'
CGG      INTEGER, PARAMETER:: NNNW = 120

      POINTER (PNTRIQ,RIQDUMMY)
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB10/NH(NNNW)
      COMMON/iounit/istdi,istdo,istde

      GOTO 2

    1 WRITE (istde,*) ' redo ...'
*
*   Read a record
*
    2 NSUBS = 0
      READ (*,'(A)') RECORD

      WRITE(734,'(a)') trim(record)

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
            IEND = I - 1
*
*   Parse the substring from left to right
*
*   (1) Determine the principal quantum number
*
            IF (RECORD(ISTART:ISTART) .EQ. '*') THEN
               NPANY = .TRUE.
               ISTART = MIN (ISTART + 1, IEND)
            ELSE
               NPANY = .FALSE.
               NFORM = 0
               DO J = ISTART, IEND
                  IF (LDIGIT (RECORD(J:J))) NFORM = NFORM + 1
               ENDDO
               IF (NFORM .EQ. 0) THEN
                  WRITE (istde,*) 'GETRSL: Unable to interpret '
     &                         , 'the principal quantum number;'
                  GOTO 1
               ENDIF
               CALL CONVRT (NFORM, CNUM, LENTH)
               FORM = '(1I'//CNUM(1:LENTH)//')'
               READ (RECORD(ISTART:ISTART+NFORM-1),
     :               FORM,IOSTAT = IOS) NPQ
               IF (IOS .NE. 0) THEN
                  WRITE (istde,*) 'GETRSL: Unable to interpret ',
     & 'string ',RECORD(ISTART:IEND-2),' as a principal quantum number'
                  GOTO 1
               ENDIF
               ISTART = ISTART + NFORM
            ENDIF
*
*   (2) Determine the orbital angular momentum quantum number
*
            NLQ = RECORD(ISTART:ISTART)
            IF (NLQ .EQ. '*') THEN
               NLANY = .TRUE.
            ELSE
               NLANY = .FALSE.
            ENDIF
*
*   (3) Determine the spin-orbit component
*
            IF (IEND .GT. ISTART) THEN
               NSQ = RECORD(IEND:IEND)
               IF (NSQ .EQ. '*') THEN
                  NSANY = .TRUE.
               ELSEIF (NSQ .EQ. '-') THEN
                  NSANY = .FALSE.
               ELSE
                  WRITE (istde,*) 'GETRSL: Unable to interpret '
     & ,'string ', NSQ, ' as a spin-orbit component indicator'
                  GOTO 1
               ENDIF
            ELSE
               IF (NLANY) THEN
                  NSANY = .TRUE.
               ELSE
                  NSANY = .FALSE.
                  NSQ = ' '
               ENDIF
            ENDIF
*
* 
            FOUND = .FALSE.
            DO J = 1, NW
               NPOK = NPANY .OR. (NP(J) .EQ. NPQ)
               NLOK = NLANY .OR. (NLQ .EQ. NH(J)(1:1))
               NSOK = NSANY .OR. (NSQ .EQ. NH(J)(2:2))
               IF (NPOK .AND. NLOK .AND. NSOK) THEN
                  DO K = 1, NSUBS
                     IF (indx(K) .EQ. J) THEN
                        WRITE (istde,*) 'GETRSL: ',
     &                                 'Repeated subshell in list;'
                        GOTO 1
                     ENDIF
                  ENDDO
                  FOUND = .TRUE.
                  NSUBS = NSUBS + 1
                  indx(NSUBS) = J
               ENDIF
            ENDDO
*
            IF (.NOT. FOUND) THEN
               WRITE (istde,*) 'GETRSL: Subshell not occupied as '
     &                      , ' according to CSL  File;'
               GOTO 1
            ENDIF
*
            ISTART = 0
         ENDIF
      ENDIF
*
      IF (I .LT. 500) THEN
         I = I + 1
         GOTO 3
      ENDIF
*
      RETURN
      END
