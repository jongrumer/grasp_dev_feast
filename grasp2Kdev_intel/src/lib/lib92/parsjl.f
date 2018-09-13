************************************************************************
*                                                                      *
      SUBROUTINE PARSJL (MODE,NCORE,RECORD,LOC,JX,NJX,IERR)
*                                                                      *
*   READs and  parses a string that specifies angular momentum quan-   *
*   tum numbers.                                                       *
*                                                                      *
*   Call(s) to: [LIB92] CONVRT.                                        *
*                                                                      *
*   Written by Farid A Parpia              Last revised: 21 Dec 1992   *
*                                                                      *
************************************************************************
*

! This subroutine is here to replace the existing one which has been 
! renamed as PARSJL_OLD. The purpose is to remove illegal GOTO into an 
! IF block. The strategy is simple: copy the kernl to the bottom
! and add a do-loop over I.
! To restore (re-use) this subroutine, give the existing PARSJL a new 
! name and then change the name of this subroutine back to PARSJL, and 
! compile.
! XHH 1997.01.28

Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 CTEGER
      CHARACTER*1 RECI
*
      DIMENSION JX(*)
*
      COMMON/ORB2/NCF,NW,PNTRIQ
*
*   There cannot be more than JXMAX angular momenta specified; if
*   MODE is 1, the subshell quantum numbers are being read; if MODE
*   is 2, the intermediate and final angular momenta are being read
*
      IF (MODE .EQ. 1) THEN
         NJXMAX = 2*(NW-NCORE)
      ELSE
         NJXMAX = NW-NCORE
      ENDIF
*
*   Initialise NJX
*
      NJX = 0
*
*   Parse RECORD from left to right
*
      ISTART = 0
      I = 1

! The original algorithm goes through the whole subroutine at least
! once, whatever the value of LOC. Thus we define another integer
! iloop to achieve this.
! XHH 1997.01.28

      iloop = MAX(1,LOC)
      DO I = 1, iloop

      RECI = RECORD(I:I)
      IF ((RECI .NE. ' ') .AND.
     :    (RECI .NE. ',') .AND.
     :    (RECI .NE. ';')) THEN
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
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            NJX = NJX+1
            IF (NJX .GT. NJXMAX) THEN
               PRINT *, 'PARSJL: Too many angular momentum'
     &, ' quantum numbers specified;'
               IERR = 1
               GOTO 3
            ENDIF
            IEND = I-1
            IF (IFRAC .EQ. 0) THEN
               CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IEND),FMT = FORM) J
               IF (J .LT. 0) THEN
                  PRINT *, 'PARSJL: Negative angular momentum'
     &, ' quantum number found;'
                  IERR = 2
                  GOTO 3
               ENDIF
               JX(NJX) = 2*J
            ELSE
               CALL CONVRT (IEND-IFRAC,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(IFRAC+1:IEND),FMT = FORM) J
               IF (J .NE. 2) THEN
                  PRINT *, 'PARSJL: The denominator of a'
     &, ' fractional quantum number must be 2;'
                  IERR = 3
                  GOTO 3
               ENDIF
               CALL CONVRT (IFRAC-ISTART,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IFRAC-1),FMT = FORM) J
               IF (J .LT. 0) THEN
                  PRINT *, 'PARSJL: Negative angular momentum'
     &, ' quantum number found;'
                  IERR = 4
                  GOTO 3
               ENDIF
               JX(NJX) = J
            ENDIF
            ISTART = 0
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         ENDIF
      ENDIF

      ENDDO

      IF (LOC .LT. 1) GOTO 10

! The following was accessed one extra time only when 1 <= I = LOC+1 .
! After the do-loop above, the value of I would be either LOC+1 or
! 2, depending on if LOC is greater than 1 or not, respectively


!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! The following is exactly the same as those above 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            NJX = NJX+1
            IF (NJX .GT. NJXMAX) THEN
               PRINT *, 'PARSJL: Too many angular momentum'
     &, ' quantum numbers specified;'
               IERR = 1
               GOTO 3
            ENDIF
            IEND = I-1
            IF (IFRAC .EQ. 0) THEN
               CALL CONVRT (IEND-ISTART+1,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IEND),FMT = FORM) J
               IF (J .LT. 0) THEN
                  PRINT *, 'PARSJL: Negative angular momentum'
     &, ' quantum number found;'
                  IERR = 2
                  GOTO 3
               ENDIF
               JX(NJX) = 2*J
            ELSE
               CALL CONVRT (IEND-IFRAC,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(IFRAC+1:IEND),FMT = FORM) J
               IF (J .NE. 2) THEN
                  PRINT *, 'PARSJL: The denominator of a'
     &, ' fractional quantum number must be 2;'
                  IERR = 3
                  GOTO 3
               ENDIF
               CALL CONVRT (IFRAC-ISTART,CTEGER,LTEGER)
               FORM = '(1I'//CTEGER(1:LTEGER)//')'
               READ (RECORD(ISTART:IFRAC-1),FMT = FORM) J
               IF (J .LT. 0) THEN
                  PRINT *, 'PARSJL: Negative angular momentum'
     &, ' quantum number found;'
                  IERR = 4
                  GOTO 3
               ENDIF
               JX(NJX) = J
            ENDIF
            ISTART = 0
!XHH~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*
   10 IERR = 0
*
    3 RETURN
      END
