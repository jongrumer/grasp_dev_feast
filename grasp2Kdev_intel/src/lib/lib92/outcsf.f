************************************************************************
*                                                                      *
      SUBROUTINE OUTCSF (ICSF,NCORE,NW,IQ,ISPAR,ITJPO,JCUP,JQS)
*                                                                      *
*   Writes CSFs to .csf file.                                          *
*                                                                      *
*   Call(s) to: >EXTERNAL<: IQ, ISPAR, ITJPO, JCUP, JQS.               *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      CHARACTER*256 RECORD
      CHARACTER*6 FORM
      CHARACTER*2 NH
      CHARACTER*1 CNUM
*
      EXTERNAL IQ,ISPAR,ITJPO,JCUP,JQS

      DIMENSION CNUM(2)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNW = 214)
      COMMON/ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
      DATA CNUM /'1','2'/
*
*
*   Write out the peel subshell information for occupied subshells
*   only
*
      NOC = 0
      IEND = 0
      DO 1 II = NCORE+1,NW
         IOCII = IQ (II,ICSF)
         IF (IOCII .GT. 0) THEN
            NOC = NOC+1
            IBEG = IEND+1
            IEND = IBEG+8
            WRITE (RECORD(IBEG:IEND),300) NP(II),NH(II),IQ (II,ICSF)
         ENDIF
    1 CONTINUE
      WRITE (21,'(A)') RECORD(1:IEND)
*
*   Write out the subshell total angular momentum only if the
*   subshell is open; write out the seniority only if j = 7/2,
*   q = 4, and J = 2 or 4
*
      IEND = 0
      DO 2 II = NCORE+1,NW
         IOCII = IQ (II,ICSF)
         IF (IOCII .GT. 0) THEN
            N2JII = NKJ(II)
            IF (IOCII .LT. N2JII+1) THEN
               JTWICE = JQS (3,II,ICSF)-1
               IF ((N2JII .EQ. 7) .AND. (IOCII .EQ. 4) .AND.
     :             ((JTWICE .EQ. 4) .OR. (JTWICE .EQ. 8))) THEN
                  IBEG = IEND+1
                  IEND = IBEG+2
                  RECORD(IBEG:IEND) = '   '
                  IBEG = IEND+1
                  IEND = IBEG
                  WRITE (RECORD(IBEG:IEND),'(1I1)')
     :               JQS (1,II,ICSF)
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
*
*   Intra-subshell coupling information
*
      RECORD(1:9) = '         '
      IEND = 9
      IOC = 0
      IOPEN = 0
      DO 3 II = NCORE+1,NW-1
         IOCII = IQ (II,ICSF)
         IF (IOCII .GT. 0) THEN
            IOC = IOC+1
            JTWICE = JQS (3,II,ICSF)-1
            IF (IOC .LT. NOC) THEN
               IF (NKJ(II)+1-IOCII .GT. 0) IOPEN = IOPEN+1
               IBEG = IEND+1
               IEND = IBEG+8
               RECORD(IBEG:IEND) = '         '
               IF ((IOPEN .GE. 2) .AND. (JTWICE .GT. 0))THEN
                  JTWICE = JCUP (IOPEN-1,ICSF)-1
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
      JTWICE = ITJPO (ICSF)-1
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
      IF (ISPAR (ICSF) .EQ. 1) THEN
         RECORD(IBEG:IEND) = '+'
      ELSE
         RECORD(IBEG:IEND) = '-'
      ENDIF
      WRITE (21,'(A)') RECORD(1:IEND)
*
      RETURN
*
  300 FORMAT (1X,1I2,1A2,'(',1I2,')')
*
      END
