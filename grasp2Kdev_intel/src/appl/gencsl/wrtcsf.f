************************************************************************
*                                                                      *
      SUBROUTINE WRTCSF (NCORE,NORB,IOCCS,JATOM,IPTY,PRNTSC)
*                                                                      *
*   Writes CSFs to grasp92.csf; echoes to screen if PRNTSC is .TRUE.   *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 16 Oct 1994   *
*                                                                      *
************************************************************************
*
      LOGICAL PRNTSC
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 SYM
      CHARACTER*1 CNUM
*
      PARAMETER (LJLMAX = 20)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWM1
CGG      PARAMETER (NNNWM1 = 119)
CGG      INTEGER NNNWM2
CGG      PARAMETER (NNNWM2 = 118)
*
      DIMENSION IOCCS(NNNW),CNUM(2)
*
*
      POINTER (PJCLST,JCLIST(NNNWM1,*))
*
      COMMON/COUBOX/PJCLST,MNJVC,LJCL(NNNWM1),ICPTR(NNNWM2)
     :      /ORBBOX/JVLIST(NNNW,LJLMAX),JWLIST(NNNW,LJLMAX),
     :              JLIST(NNNW,LJLMAX),LJL(NNNW),IOPTR(NNNW)
     :      /ORBNUM/NP(NNNW),N2J(NNNW),NL(NNNW)
     :      /ORBSYM/SYM(NNNW)
*
      DATA CNUM /'1','2'/
*
*
*   Write out the peel subshell information for occupied subshells
*   only
*
      NOC = 0
      IEND = 0
      DO 1 II = NCORE+1,NORB
         IOCII = IOCCS(II)
         IF (IOCII .GT. 0) THEN
            NOC = NOC+1
            IBEG = IEND+1
            IEND = IBEG+8
            WRITE (RECORD(IBEG:IEND),300) NP(II),SYM(II),IOCCS(II)
         ENDIF
    1 CONTINUE
      WRITE (21,'(A)') RECORD(1:IEND)
      IF (PRNTSC) WRITE (*,'(A)') RECORD(1:IEND)
*
*   Write out the subshell total angular momentum only if the
*   subshell is open; write out the seniority only if j = 7/2,
*   q = 4, and J = 2 or 4
*
      IEND = 0
      DO 2 II = NCORE+1,NORB
         IOCII = IOCCS(II)
         IF (IOCII .GT. 0) THEN
            N2JII = N2J(II)
            IF (IOCII .LT. N2JII+1) THEN
               JTWICE = JLIST(II,IOPTR(II))
               IF ((N2JII .EQ. 7) .AND. (IOCII .EQ. 4) .AND.
     :             ((JTWICE .EQ. 4) .OR. (JTWICE .EQ. 8))) THEN
                  IBEG = IEND+1
                  IEND = IBEG+2
                  RECORD(IBEG:IEND) = '   '
                  IBEG = IEND+1
                  IEND = IBEG
                  WRITE (RECORD(IBEG:IEND),'(1I1)')
     :               JVLIST(II,IOPTR(II))
                  IBEG = IEND+1
                  IEND = IBEG
                  RECORD(IBEG:IEND) = ';'
               ELSE
                  IBEG = IEND+1
                  IEND = IBEG+4
                  RECORD(IBEG:IEND) = '     '
               ENDIF
               IBEG = IEND+1
               IF (MOD (JTWICE,2) .EQ. 0) THEN
                  IEND = IBEG+3
                  WRITE (RECORD(IBEG:IEND),'(2X,1I2)') JTWICE/2
               ELSE
                  IEND = IBEG+1
                  WRITE (RECORD(IBEG:IEND),'(1I2)') JTWICE
                  IBEG = IEND+1
                  IEND = IBEG+1
                  RECORD(IBEG:IEND) = '/2'
               ENDIF
            ELSE
               IBEG = IEND+1
               IEND = IBEG+8
               RECORD(IBEG:IEND) = '         '
            ENDIF
         ENDIF
    2 CONTINUE
      WRITE (21,'(A)') RECORD(1:IEND)
      IF (PRNTSC) WRITE (*,'(A)') RECORD(1:IEND)
*
*   Intra-subshell coupling information
*
      RECORD(1:9) = '         '
      IEND = 9
      IOC = 0
      IOPEN = 0
      DO 3 II = NCORE+1,NORB-1
         IOCII = IOCCS(II)
         IF (IOCII .GT. 0) THEN
            IOC = IOC+1
            JTWICE = JLIST(II,IOPTR(II))
            IF (IOC .LT. NOC) THEN
               IF (N2J(II)+1-IOCII .GT. 0) IOPEN = IOPEN+1
               IBEG = IEND+1
               IEND = IBEG+8
               RECORD(IBEG:IEND) = '         '
               IF ((IOPEN .GE. 2) .AND. (JTWICE .GT. 0))THEN
                  ICLOC = II-NCORE-1
                  JTWICE = JCLIST(ICLOC,ICPTR(ICLOC))
                  IF (JTWICE .NE. 0) THEN
                     IF (MOD (JTWICE,2) .EQ. 0) THEN
                        LENTH = INT (LOG10 (REAL (JTWICE/2)))+1
                        FORM = '(1I'//CNUM(LENTH)//')'
                        WRITE (RECORD(IBEG:IBEG+LENTH-1),FORM) JTWICE/2
                     ELSE
                        LENTH = INT (LOG10 (REAL (JTWICE)))+1
                        FORM = '(1I'//CNUM(LENTH)//')'
                        WRITE (RECORD(IBEG:IBEG+LENTH-1),FORM) JTWICE
                        IBEG = IBEG+LENTH
                        RECORD(IBEG:IBEG+1) = '/2'
                     ENDIF
                  ELSE
                     RECORD(IBEG:IBEG+1) = '0'
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
    3 CONTINUE
*
*   Total angular momentum and parity
*
      IBEG = IEND+1
      IEND = IBEG+4
      RECORD(IBEG:IEND) = '     '
      JTWICE = JATOM
      IF (JTWICE .NE. 0) THEN
         IF (MOD (JTWICE,2) .EQ. 0) THEN
            LENTH = INT (LOG10 (REAL (JTWICE/2)))+1
            FORM = '(1I'//CNUM(LENTH)//')'
            WRITE (RECORD(IBEG:IBEG+LENTH-1),FORM) JTWICE/2
         ELSE
            LENTH = INT (LOG10 (REAL (JTWICE)))+1
            FORM = '(1I'//CNUM(LENTH)//')'
            WRITE (RECORD(IBEG:IBEG+LENTH-1),FORM) JTWICE
            IBEG = IBEG+LENTH
            RECORD(IBEG:IBEG+1) = '/2'
         ENDIF
      ELSE
         RECORD(IBEG:IBEG+1) = '0'
      ENDIF
      DO 4 I = IEND-4,IEND
         IF (RECORD(I:I) .EQ. ' ') THEN
            IBEG = I
            GOTO 5
         ENDIF
    4 CONTINUE
    5 IEND = IBEG
      IF (IPTY .EQ. 1) THEN
         RECORD(IBEG:IEND) = '+'
      ELSE
         RECORD(IBEG:IEND) = '-'
      ENDIF
      WRITE (21,'(A)') RECORD(1:IEND)
      IF (PRNTSC) WRITE (*,'(A)') RECORD(1:IEND)
*
      RETURN
*
  300 FORMAT (1X,1I2,1A2,'(',1I2,')')
*
      END
