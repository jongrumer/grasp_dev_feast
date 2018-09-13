************************************************************************
*                                                                      *
      SUBROUTINE LODCSLM (NCORE)
*                                                                      *
*   Loads the data from the  .csl  file. A number of checks are made   *
*   to ensure correctness and consistency.                             *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
*                        PACK, PARSJL, PRSRCN, PRSRSL, RALC2D,         *
*                        RALLOC.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*   Updated by Xinghong He                               31 Oct 1997
*       To accept both block and non-block formats
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H,O-Z)
      include 'parameters.def'
CFF      PARAMETER (NNNW = 120)
CFF      PARAMETER (NNNWP = 30)
      LOGICAL EMPTY,FULL,LDBPA
      CHARACTER*256 RECORD
      CHARACTER*2 NH
      CHARACTER*1 RECL
*
      PARAMETER (NW2 = 2*NNNW)
*
      DIMENSION IOCC(NNNW),IQSUB(NW2),JX(NNNW)
*
      INTEGER*4  IQA, JQSA, JCUPA
      POINTER (PNTRIQ,IQA(NNNWP,1))
      POINTER (PNTJQS,JQSA(NNNWP,3,1))
      POINTER (PNJCUP,JCUPA(NNNWP,1))
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
*
*   Entry message
*
      PRINT *, 'Loading Configuration Symmetry List File ...'
*
*   Get the list of subshells
*
      NW = 0
*
*   Read the list of core subshells; set up the arrays NP, NAK,
*   NKL, NKJ, NH for these subshells
*
      CALL PRSRSL (21,1)
      NCORE = NW
      NCORP1 = NW+1
*
*   Skip the peel subshell identification header; read the list of
*   peel subshells; set up the arrays NP, NAK, NKL, NKJ, NH for
*   these subshells
*
      READ (21,*)
      CALL PRSRSL (21,2)
      NPEEL = NW-NCORE
      NPEEL2 = NPEEL*2
*
*   Ensure that the sets of core and peel subshell are disjoint
*
      DO 2 J = NCORE+1,NW
         NPJ = NP(J)
         NAKJ = NAK(J)
         DO 1 I = 1,NCORE
            IF ((NP(I) .EQ. NPJ) .AND. (NAK(I) .EQ. NAKJ)) THEN
               PRINT *, 'LODCSL: The lists of core and'
               PRINT *, ' peel subshells must form'
               PRINT *, ' disjoint sets.'
               STOP
            ENDIF
    1    CONTINUE
    2 CONTINUE
*
*   Print the number of relativistic subshells
*
      IF (NW .GT. 1) THEN
         CALL CONVRT (NW,RECORD,LENTH)
         PRINT *, ' there are '//RECORD(1:LENTH)
     :          //' relativistic subshells;'
      ELSE
         PRINT *, ' there is 1 relativistic subshell;'
      ENDIF
*
*   Initial allocation for arrays with a dimension dependent
*   on the number of CSFs; the initial allocation must be
*   greater than 1
*
      NCFD = 2
      CALL ALLOC (PNTRIQ,NNNWP  *NCFD,4)
      CALL ALLOC (PNTJQS,NNNWP*3*NCFD,4)
      CALL ALLOC (PNJCUP,NNNWP  *NCFD,4)
*
*   Skip the header for the list of CSFs
*
      READ (21,*)
*
*   NREC is the sequence number of the last record read in the
*   Configuration Symmetry List File
*
      NREC = 5
*
*   There must be three records for each CSF: For instance,
*
*    4s ( 2) 4p-( 2) 4p ( 3) 4d-( 2) 4d ( 5) 4f-( 6) 4f ( 4)
*                        3/2       0     5/2             2;4
*                                           1               3-
*   Each record is as follows:
*      (1) Peel subshell occupation number (q) specification.
*      (2) Subshell total angular momentum quantum number
*          (v_, J_sub) specifications. Subshell total angular
*          momentum quantum numbers (J_sub) must be specified
*          for all open subshells, even if this can be
*          deduced from the subshell occupation. Seniority
*          quantum numbers (v) must be specified for
*          subshells with j = 7/2, q = 4, J_sub = 2 or 4.
*          The seniority quantum number must precede the
*          total angular momentum quantum number: in the
*          example above, for the 4f_7/2 subshell, q = 4
*          and J_sub = 4, whence it is necessary to specify
*          v --- 2 in this case.
*      (3) The minimum number of intermediate angular
*          momentum quantum numbers (X) as well as the final
*          angular momentum quantum number (J) immediately
*          followed by the sign of the parity (P) must be
*          specified on this record. In the example above,
*          the J_sub = 0 for the 4d_3/2 subshell, whence
*          it is unnecessary to specify its coupling to all
*          preceding subshells.
*   These conventions have been chosen so as to render the CSF
*   specifications easily interpreted by the user
*
      NCF = 0
    3 NCF = NCF+1
*
      READ (21,'(A)',IOSTAT = IOS) RECORD
***********************************************************************
*blk*
*   To skip the border line added to mark the end of a block
*
      IF (RECORD(1:2) .EQ. ' *') THEN
         READ (21,'(A)',IOSTAT = IOS) RECORD
      ENDIF
***********************************************************************

      IF (IOS .EQ. 0) THEN
*
*   Read in the occupations (q) of the peel shells; stop with a
*   message if an error occurs
*
         CALL PRSRCN (RECORD,NCORE,IOCC,IERR)
         IF (IERR .NE. 0) GOTO 26
*
*   Read the J_sub and v quantum numbers
*
         READ (21,'(A)',IOSTAT = IOS) RECORD
         IF (IOS .NE. 0) THEN
            PRINT *, 'LODCSL: Expecting subshell quantum'
            PRINT *, ' number specification;'
            GOTO 26
         ENDIF
         LOC = len_trim (RECORD)
         CALL PARSJL (1,NCORE,RECORD,LOC,IQSUB,NQS,IERR)
         IF (IERR .NE. 0) GOTO 26
*
*   Read the X, J, and (sign of) P quantum numbers
*
         READ (21,'(A)',IOSTAT = IOS) RECORD
         IF (IOS .NE. 0) THEN
            PRINT *, 'LODCSL: Expecting intermediate'
            PRINT *, ' and final angular momentum'
            PRINT *, ' quantum number and final parity'
            PRINT *, ' specification;'
            GOTO 26
         ENDIF
*
*   Allocate additional storage if necessary
*
         IF (NCF .GT. NCFD) THEN
            NEWSIZ = NCFD+NCFD/2
            CALL RALC2D (PNTRIQ,NNNWP  ,NCFD,NNNWP  ,NEWSIZ,4)
            CALL RALC2D (PNTJQS,NNNWP*3,NCFD,NNNWP*3,NEWSIZ,4)
            CALL RALC2D (PNJCUP,NNNWP  ,NCFD,NNNWP  ,NEWSIZ,4)
            NCFD = NEWSIZ
         ENDIF
*
*   Zero out the arrays that store packed integers
*
         DO 4 I = 1,NNNWP
            IQA(I,NCF) = 0
            JQSA(I,1,NCF) = 0
            JQSA(I,2,NCF) = 0
            JQSA(I,3,NCF) = 0
            JCUPA(I,NCF) = 0
    4    CONTINUE
*
*   Determine the parity and all intermediate and the final
*   angular momentum quantum numbers
*
         DO 5 I = 256,1,-1
            IF (RECORD(I:I) .NE. ' ') THEN
               LOC = I
               GOTO 6
            ENDIF
    5    CONTINUE
    6    RECL = RECORD(LOC:LOC)
         IF     (RECL .EQ. '+') THEN
            ISPARC = +1
         ELSEIF (RECL .EQ. '-') THEN
            ISPARC = -1
         ELSE
            PRINT *, 'LODCSL: Incorrect parity'
            PRINT *, ' specification;'
            GOTO 26
         ENDIF
         LOC = LOC-1
*
         CALL PARSJL (2,NCORE,RECORD,LOC,JX,NJX,IERR)
         IF (IERR .NE. 0) GOTO 26
*
*   Set the occupation and subshell quantum number array elements
*   in IQ, JQS for the core subshells
*
         DO 7 I = 1,NCORE
            CALL PACK (NKJ(I)+1,I,IQA(1,NCF))
            CALL PACK (0,I,JQSA(1,1,NCF))
            CALL PACK (0,I,JQSA(1,2,NCF))
            CALL PACK (1,I,JQSA(1,3,NCF))
    7    CONTINUE
*
*   Check all subshell, intermediate and final angular momentum
*   quantum numbers; set the array elements in IQ, JQS for the peel
*   subshells; set the coupling array element in JCUP and the total
*   angular momentum array element in ITJPO
*
         IOC = 0
         IPTY = 0
         NQSN = 0
         NJXN = 0
         NPEELN = 0
         NOPEN = 0
         JLAST = 0
         ILAST = 0
         DO 12 I = NCORP1,NW
            IOCCI = IOCC(I)
            NPEELN = NPEELN+IOCCI
            NKJI = NKJ(I)
            IFULLI = NKJI+1
            EMPTY = IOCCI .EQ. 0
            IF (.NOT. EMPTY) IOC = IOC+1
            FULL = IOCCI .EQ. IFULLI
            IF (EMPTY .OR. FULL) THEN
               NU = 0
               JSUB = 0
            ELSE
               IPTY = IPTY+NKL(I)*IOCCI
               IF (NKJI .NE. 7) THEN
                  NQSN = NQSN+1
                  IF (NQSN .GT. NQS) THEN
                     PRINT *, 'LODCSL: Too few subshell'
                     PRINT *, ' quantum numbers specified;'
                     GOTO 26
                  ENDIF
                  NU = 0
                  JSUB = IQSUB(NQSN)
               ELSE
                  IF (IOCCI .NE. 4) THEN
                     NQSN = NQSN+1
                     IF (NQSN .GT. NQS) THEN
                        PRINT *, 'LODCSL: Too few subshell'
                        PRINT *, ' quantum numbers specified;'
                        GOTO 26
                     ENDIF
                     NU = 0
                     JSUB = IQSUB(NQSN)
                  ELSE
                     NQSN = NQSN+1
                     IF (NQSN .GT. NQS) THEN
                        PRINT *, 'LODCSL: Too few subshell'
                        PRINT *, ' quantum numbers specified;'
                        GOTO 26
                     ENDIF
                     JSUB = IQSUB(NQSN)
                     IF ((JSUB .EQ. 4) .OR. (JSUB .EQ. 8)) THEN
                        NU = JSUB/2
                        NQSN = NQSN+1
                        IF (NQSN .GT. NQS) THEN
                           PRINT *, 'LODCSL: Too few subshell'
                           PRINT *, ' quantum numbers specified;'
                           GOTO 26
                        ENDIF
                        JSUB = IQSUB(NQSN)
                     ELSE
                        NU = 0
                     ENDIF
                  ENDIF
               ENDIF
               IQT = MIN (IOCCI,IFULLI-IOCCI)
               LOC = (IFULLI-2)/2
               LOC = (LOC*(LOC+1))/2+IQT
               NBEG = JTAB(LOC+1)+1
               NEND = JTAB(LOC+2)
               DO 8 J = NBEG,NEND,3
                  IF (NTAB(J+2) .EQ. JSUB+1) THEN
                     IF (NU .EQ. 0) THEN
                        NU = NTAB(J)
                        GOTO 9
                     ELSE
                        IF (NTAB(J) .EQ. NU) GOTO 9
                     ENDIF
                  ENDIF
    8          CONTINUE
               CALL CONVRT (NP(I),RECORD,LENTH)
               PRINT *, 'LODCSL: Subshell quantum numbers'
               PRINT *, ' specified incorrectly for'
               PRINT *, ' '//RECORD(1:LENTH)//NH(I)//' subshell.'
               GOTO 26
            ENDIF
    9       IF ((.NOT. EMPTY) .AND. (.NOT. FULL)) THEN
               NOPEN = NOPEN+1
               IF (NOPEN .GT. 1) THEN
                  IF (JSUB .EQ. 0) THEN
                     JXN = JLAST
                  ELSE
                     ILAST = IOC
                     NJXN = NJXN+1
                     IF (NJXN .GT. NJX) THEN
                        PRINT *, 'LODCSL: Too few intermediate'
                        PRINT *, ' and final angular momentum'
                        PRINT *, ' quantum numbers specified;'
                        GOTO 26
                     ENDIF
                     JXN = JX(NJXN)
                     DO 10 J = ABS (JLAST-JSUB),JLAST+JSUB,2
                        IF (JXN .EQ. J) GOTO 11
   10                CONTINUE
                     CALL CONVRT (NP(I),RECORD,LENTH)
                     PRINT *, 'LODCSL: coupling of '
     :                      //RECORD(1:LENTH)//NH(I)
                     PRINT *, ' subshell to previous subshells'
                     PRINT *, ' is incorrect.'
                     GOTO 26
                  ENDIF
   11             CALL PACK (JXN+1,NOPEN-1,JCUPA(1,NCF))
                  JLAST = JXN
               ELSE
                  JLAST = JSUB
               ENDIF
            ENDIF
            CALL PACK (IOCCI,I,IQA(1,NCF))
            CALL PACK (NU    ,I,JQSA(1,1,NCF))
            CALL PACK (0     ,I,JQSA(1,2,NCF))
            CALL PACK (JSUB+1,I,JQSA(1,3,NCF))
   12    CONTINUE
*
         DO 13 I = MAX (1,NOPEN),NW
            CALL PACK (0,I,JCUPA(1,NCF))
   13    CONTINUE
*
         IF (NQSN .NE. NQS) THEN
            PRINT *, 'LODCSL: Too many subshell'
            PRINT *, ' quantum numbers specified;'
            GOTO 26
         ENDIF
*
         IF (ILAST .NE. IOC) NJXN = NJXN+1
         IF (NJXN .NE. NJX) THEN
            PRINT *, 'LODCSL: Too many intermediate'
            PRINT *, ' and final angular momentum'
            PRINT *, ' quantum numbers specified;'
            GOTO 26
         ENDIF
*
         IF (JX(NJXN) .NE. JLAST) THEN
            PRINT *, 'LODCSL: Final angular momentum'
            PRINT *, ' incorrectly specified;'
            GOTO 26
         ENDIF
*
         IPTY = (-1)**IPTY
         IF (IPTY .NE. ISPARC) THEN
            PRINT *, 'LODCSL: Parity specified incorrectly;'
            GOTO 26
         ENDIF
*
         JPI = (JLAST+1)*IPTY
         CALL PACK (JPI,NNNW,JCUPA(1,NCF))
*
         IF (NCF .GT. 1) THEN
            IF (NPEELN .NE. NPEEL) THEN
               PRINT *, 'LODCSL: Inconsistency in the number'
               PRINT *, ' of electrons.'
               GOTO 26
            ENDIF
         ELSE
            NPEEL = NPEELN
         ENDIF
*
*   Check if this CSF was already in the list; stop with a
*   message if this is the case
*
C         print *, 'Check duplicated CSFs'
         IF (NCF .GT. 1) THEN
            DO 16 J = 1,NCF-1
C                  print *,'j= ',j,ncf
               DO 14 I = NCORP1,NW
C                  print *,i
C                  print *, IQ(I,J), JQS(1,I,J), JQS(2,I,J),JQS(3,I,J)
C          print *, IQ(I,ncf), JQS(1,I,ncf), JQS(2,I,ncf),JQS(3,I,ncf)
                  IF (IQ (I,J) .NE. IQ (I,NCF)) GOTO 17
                  IF (JQS (1,I,J) .NE. JQS (1,I,NCF)) GOTO 17
                  IF (JQS (2,I,J) .NE. JQS (2,I,NCF)) GOTO 17
                  IF (JQS (3,I,J) .NE. JQS (3,I,NCF)) GOTO 17
   14          CONTINUE
               DO 15 I = 1,NOPEN-1
                  print *,i
                  print *, JCUP (I,J),JCUP (I,NCF)
                  IF (JCUP (I,J) .NE. JCUP (I,NCF)) GOTO 17
   15          CONTINUE
   16       CONTINUE
            PRINT *, 'LODCSL: Repeated CSF;'
            GOTO 26
         ENDIF
*
*   Successfully read a CSF; update NREC and read another CSF
*
   17    NREC = NREC+3
         GOTO 3
*
      ELSE
*
*   There is always at least one CSF
*
         IF (NCF .EQ. 1) THEN
            DO 18 I = 1,NCORE
               CALL PACK (NKJ(I)+1,I,IQA(1,1))
               CALL PACK (0,I,JQSA(1,1,1))
               CALL PACK (0,I,JQSA(1,2,1))
               CALL PACK (1,I,JQSA(1,3,1))
   18       CONTINUE
            CALL PACK (0,1,JCUPA(1,1))
            CALL PACK (1,NNNW,JCUPA(1,1))
         ELSE
            NCF = NCF-1
         ENDIF
*
      ENDIF
*
*   Check if any subshell is empty; eliminate it from the
*   list if this is the case; issue a message
*
      I = NCORP1
   19 IF (I .LE. NW) THEN
         DO 20 J = 1,NCF
            IF (IQ (I,J) .NE. 0) GOTO 23
   20    CONTINUE
         CALL CONVRT (NP(I),RECORD,LENTH)
         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is empty'
         PRINT *, ' in all CSFs; eliminating this'
         PRINT *, ' subshell from the list;'
         NW = NW-1
         DO 22 II = I,NW
            NP(II) = NP(II+1)
            NAK(II) = NAK(II+1)
            NKL(II) = NKL(II+1)
            NKJ(II) = NKJ(II+1)
            NH(II) = NH(II+1)
            DO 21 J = 1,NCF
               ITEMP = IQ (II+1,J)
               CALL PACK (ITEMP,II,IQA(1,J))
               ITEMP = JQS (1,II+1,J)
               CALL PACK (ITEMP,II,JQSA(II,1,J))
               ITEMP = JQS (2,II+1,J)
               CALL PACK (ITEMP,II,JQSA(II,2,J))
               ITEMP = JQS (3,II+1,J)
               CALL PACK (ITEMP,II,JQSA(II,3,J))
   21       CONTINUE
   22    CONTINUE
   23    I = I+1
         GOTO 19
      ENDIF
*
*   Store the number of electrons in the COMMON variable
*
      NCOREL = 0
      DO 24 I = 1,NCORE
         NCOREL = NCOREL+NKJ(I)+1
   24 CONTINUE
      NELEC = NCOREL+NPEEL
*
*   All done; report
*
      CALL CONVRT (NCF,RECORD,LENTH)
      PRINT *, ' there are '//RECORD(1:LENTH)
     :         //' relativistic CSFs;'
      PRINT *, ' ... load complete;'
*
*   Debug printout
*
      IF (LDBPA(1)) THEN
         WRITE (99,*) 'From LODCSL:'
         DO 25 I = 1,NCF
            WRITE (99,*) 'CSF ',I
            WRITE (99,*) 'ITJPO: ',ITJPO (I)
            WRITE (99,*) 'ISPAR: ',ISPAR (I)
            WRITE (99,*) 'IQ: ',(IQ (J,I),J = 1,NW)
            WRITE (99,*) 'JQS(1): ',(JQS (1,J,I),J = 1,NW)
            WRITE (99,*) 'JQS(2): ',(JQS (2,J,I),J = 1,NW)
            WRITE (99,*) 'JQS(3): ',(JQS (3,J,I),J = 1,NW)
            WRITE (99,*) 'JCUP: ',(JCUP (J,I),J = 1,NW-1)
   25    CONTINUE
      ENDIF
*
      RETURN
*
   26 CALL CONVRT (NCF,RECORD,LENTH)
      PRINT *, ' CSF sequence number: '//RECORD(1:LENTH)//':'
      REWIND (21)
      DO 27 I = 1,NREC
         READ (21,*)
   27 CONTINUE
      DO 28 I = 1,3
         READ (21,'(A)',ERR = 29,END = 29) RECORD
         LENTH = len_trim (RECORD)
         PRINT *, RECORD(1:LENTH)
   28 CONTINUE
   29 CLOSE (21)
      STOP
*
      END
