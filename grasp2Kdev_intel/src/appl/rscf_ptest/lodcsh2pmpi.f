************************************************************************
*                                                                      *
      SUBROUTINE LODCSH2pmpi (nfile, ncore, jb)
      IMPLICIT REAL*8          (A-H,O-Z)

      INTEGER, PARAMETER:: LOADALL = -119    ! arbitrary non-positive
*
* IMPORTANT:
* ==========
*   If jb == LOADALL, then it loads ALL blocks. In this case, the 
*   NCFblock in the common should be the total ncf and allocations
*   should be changed to this ncf. Both should be done before calling
*   this routine - like the genuine block case.
*
*   Loads the data from the  .csl  file for the current block          *
*   All parameters are inputs:
*
*     nfile:    The unit number, usually 21
*     jb:       The block number (see above)
*     ncore:    The number of the core (sub-)shells.
*   Since NCFblock is known in the block version, allocation is done
*   once outside this routine. NCFblock is checked against the one (NCF)
*   obtained here.
*
*   Call(s) to: [LIB92]: CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,
*                        PACK, PARSJL, PRSRCN, PRSRSL
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*   Modified by C. F. Fischer for block calculation      22 May 1997   *
*   Updated by Xinghong He                               08 Jul 1998   *
************************************************************************
*
      CHARACTER(LEN=*), PARAMETER:: MYNAME='LODCSH2'

      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWP = 30)
      LOGICAL EMPTY,FULL,LDBPA
      CHARACTER*256 str
      CHARACTER*2 NH
      CHARACTER*1 RECL
*
      PARAMETER (NW2 = 2*NNNW)
*
      DIMENSION IOCC(NNNW),IQSUB(NW2),JX(NNNW)
*
      INTEGER*4 IQA,JQSA,JCUPA
      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEF1/EMN,IONCTY,NELEC,Z
     :      /ORB2/NCFblock,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP
     :      /TERMS/NROWS,ITAB(31),JTAB(32),NTAB(327)
      COMMON /mpi/ myid, nprocs, ierr
      COMMON/iounit/istdi,istdo,istde
!-----------------------------------------------------------------------

      if(myid==0) then
      IF (jb .NE. LOADALL) THEN
         PRINT *, 'Loading CSF file for block ', jb
      ELSE
         PRINT *, 'Loading CSF File for ALL blocks '
      ENDIF
      endif

      NCORP1 = NCORE + 1
      NPEEL  = NW - NCORE
* 
* NPEEL is used as 1) number of peel orbitals (here) and
*                  2) number of peel electrons (later in this routine)
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
*
*   Zero out the arrays that store packed integers - only when ncfblock>0
*
      DO jcf = 1, ncfblock
         DO I = 1, NNNWP
            IQA(I,   jcf) = 0
            JQSA(I,1,jcf) = 0
            JQSA(I,2,jcf) = 0
            JQSA(I,3,jcf) = 0
            JCUPA(I, jcf) = 0
         ENDDO
      ENDDO

      NCF = 0
    3 NCF = NCF + 1
*
      READ (nfile, '(A)', IOSTAT = IOS) str

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This IF...READ makes the routine load the entire file (all blocks)
! by ignoring the end-of-block mark

      IF (IOS .EQ. 0 .AND. str(1:2) .EQ. ' *' .AND. jb .EQ. LOADALL)
     &   READ (nfile, '(A)', IOSTAT = IOS) str
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (IOS .EQ. 0 .AND. str(1:2) .NE. ' *') THEN
*
*   Read in the occupations (q) of the peel shells; stop with a
*   message if an error occurs
*
         CALL PRSRCN (str, NCORE, IOCC, IERR)
         IF (IERR .NE. 0) GOTO 28
*
*   Read the J_sub and v quantum numbers
*
         READ (nfile,'(A)',IOSTAT = IOS) str
         IF (IOS .NE. 0) THEN
            WRITE (istde,*) MYNAME//': Expecting subshell quantum'
     &,                    ' number specification;'
            GOTO 27
         ENDIF
         LOC =  LEN_TRIM (str)
         CALL PARSJL (1, NCORE, str, LOC, IQSUB, NQS, IERR)
         IF (IERR .NE. 0) GOTO 27
*
*   Read the X, J, and (sign of) P quantum numbers
*
         READ (nfile,'(A)',IOSTAT = IOS) str
         IF (IOS .NE. 0) THEN
            WRITE (istde,*) MYNAME//': Expecting intermediate '
     &,                    'and final angular momentum'
            WRITE (istde,*) 'quantum number and final parity '
     &,                    'specification;'
            GOTO 26
         ENDIF
*
*   Zero out the arrays that store packed integers
*
         DO I = 1, NNNWP
            IQA(  I,  NCF) = 0
            JQSA( I,1,NCF) = 0
            JQSA( I,2,NCF) = 0
            JQSA( I,3,NCF) = 0
            JCUPA(I,  NCF) = 0
         ENDDO
*
*   Determine the parity and all intermediate and the final
*   angular momentum quantum numbers
*
			LOC  = LEN_TRIM (str)
    6    RECL = str(LOC:LOC)
         IF     (RECL .EQ. '+') THEN
            ISPARC = +1
         ELSEIF (RECL .EQ. '-') THEN
            ISPARC = -1
         ELSE
            WRITE (istde,*) MYNAME//': Incorrect parity '
     &,                    'specification;'
            GOTO 26
         ENDIF
         LOC = LOC - 1
*
         CALL PARSJL (2, NCORE, str, LOC, JX, NJX, IERR)
         IF (IERR .NE. 0) GOTO 26
*
*   Set the occupation and subshell quantum number array elements
*   in IQ, JQS for the core subshells
*
         DO I = 1, NCORE
            CALL PACK (NKJ(I)+1, I, IQA(1,NCF))
            CALL PACK (0,        I, JQSA(1,1,NCF))
            CALL PACK (0,        I, JQSA(1,2,NCF))
            CALL PACK (1,        I, JQSA(1,3,NCF))
         ENDDO
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
         DO 12 I = NCORP1, NW
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
                     WRITE (istde,*) MYNAME//': Too few subshell '
     &,                             'quantum numbers specified;'
                     GOTO 26
                  ENDIF
                  NU = 0
                  JSUB = IQSUB(NQSN)
               ELSE
                  IF (IOCCI .NE. 4) THEN
                     NQSN = NQSN+1
                     IF (NQSN .GT. NQS) THEN
                        WRITE (istde,*) MYNAME//': Too few subshell '
     &,                                'quantum numbers specified;'
                        GOTO 26
                     ENDIF
                     NU = 0
                     JSUB = IQSUB(NQSN)
                  ELSE
                     NQSN = NQSN+1
                     IF (NQSN .GT. NQS) THEN
                        WRITE (istde,*) MYNAME//': Too few subshell '
     &,                                'quantum numbers specified;'
                        GOTO 26
                     ENDIF
                     JSUB = IQSUB(NQSN)
                     IF ((JSUB .EQ. 4) .OR. (JSUB .EQ. 8)) THEN
                        NU = JSUB/2
                        NQSN = NQSN+1
                        IF (NQSN .GT. NQS) THEN
                           WRITE (istde,*) MYNAME//': Too few subshell'
     &,                                ' quantum numbers specified;'
                           GOTO 26
                        ENDIF
                        JSUB = IQSUB(NQSN)
                     ELSE
                        NU = 0
                     ENDIF
                  ENDIF
               ENDIF
               IQT = MIN (IOCCI, IFULLI-IOCCI)
               LOC = (IFULLI-2)/2
               LOC = (LOC*(LOC+1)) / 2 + IQT
               NBEG = JTAB(LOC+1) + 1
               NEND = JTAB(LOC+2)
               DO 8 J = NBEG, NEND, 3
                  IF (NTAB(J+2) .EQ. JSUB+1) THEN
                     IF (NU .EQ. 0) THEN
                        NU = NTAB(J)
                        GOTO 9
                     ELSE
                        IF (NTAB(J) .EQ. NU) GOTO 9
                     ENDIF
                  ENDIF
    8          CONTINUE
               CALL CONVRT (NP(I), str, LENTH)
               WRITE (istde,*) MYNAME//': Subshell quantum numbers '
     &,                       'specified incorrectly for '
     &                //str(1:LENTH)//NH(I)//' subshell.'
               GOTO 26
            ENDIF
    9       IF ((.NOT. EMPTY) .AND. (.NOT. FULL)) THEN
               NOPEN = NOPEN + 1
               IF (NOPEN .GT. 1) THEN
                  IF (JSUB .EQ. 0) THEN
                     JXN = JLAST
                  ELSE
                     ILAST = IOC
                     NJXN = NJXN + 1
                     IF (NJXN .GT. NJX) THEN
                        WRITE (istde,*) MYNAME//': Too few intermediate'
     &,                                ' and final angular momentum'
     &,                                ' quantum numbers specified;'
                        GOTO 26
                     ENDIF
                     JXN = JX(NJXN)
                     DO 10 J = ABS (JLAST-JSUB), JLAST+JSUB, 2
                        IF (JXN .EQ. J) GOTO 11
   10                CONTINUE
                     CALL CONVRT (NP(I),str,LENTH)
                     WRITE (istde,*) MYNAME//': coupling of '
     :                      //str(1:LENTH)//NH(I)
     &,                ' subshell to previous subshells is incorrect.'
                     GOTO 26
                  ENDIF
   11             CALL PACK (JXN+1, NOPEN-1, JCUPA(1,NCF))
                  JLAST = JXN
               ELSE
                  JLAST = JSUB
               ENDIF
            ENDIF
            CALL PACK (IOCCI,  I, IQA(1,NCF))
            CALL PACK (NU    , I, JQSA(1,1,NCF))
            CALL PACK (0     , I, JQSA(1,2,NCF))
            CALL PACK (JSUB+1, I, JQSA(1,3,NCF))
   12    CONTINUE
*
         DO I = MAX (1,NOPEN), NW
            CALL PACK (0, I, JCUPA(1,NCF))
         ENDDO
*
         IF (NQSN .NE. NQS) THEN
            WRITE (istde,*) MYNAME//': Too many subshell'
     &,                    ' quantum numbers specified;'
            GOTO 26
         ENDIF
*
         IF (ILAST .NE. IOC) NJXN = NJXN+1
         IF (NJXN .NE. NJX) THEN
            WRITE (istde,*) MYNAME//': Too many intermediate'
     &,                    ' and final angular momentum'
     &,                    ' quantum numbers specified;'
            GOTO 26
         ENDIF
*
         IF (JX(NJXN) .NE. JLAST) THEN
            WRITE (istde,*) MYNAME//': Final angular momentum'
     &,                    ' incorrectly specified;'
            GOTO 26
         ENDIF
*
         IPTY = (-1)**IPTY
         IF (IPTY .NE. ISPARC) THEN
            WRITE (istde,*) MYNAME//': Parity specified incorrectly;'
            GOTO 26
         ENDIF
*
         JPI = (JLAST+1)*IPTY
         CALL PACK (JPI, NNNW, JCUPA(1,NCF))
*
         IF (NCF .GT. 1) THEN
            IF (NPEELN .NE. NPEEL) THEN
               WRITE (istde,*) MYNAME//': Inconsistency in the number'
     &,                       ' of electrons.'
               GOTO 26
            ENDIF
         ELSE
            NPEEL = NPEELN
         ENDIF
*
*   Check if this CSF was already in the list; stop with a
*   message if this is the case
*
         IF (NCF .GT. 1) THEN
            DO J = 1,NCF-1
               DO I = NCORP1,NW
                  IF (IQ (I,J) .NE. IQ (I,NCF)) GOTO 17
                  IF (JQS (1,I,J) .NE. JQS (1,I,NCF)) GOTO 17
                  IF (JQS (2,I,J) .NE. JQS (2,I,NCF)) GOTO 17
                  IF (JQS (3,I,J) .NE. JQS (3,I,NCF)) GOTO 17
               ENDDO
               DO I = 1,NOPEN-1
                  IF (JCUP (I,J) .NE. JCUP (I,NCF)) GOTO 17
               ENDDO
            ENDDO
            WRITE(istde,*) MYNAME//': Repeated CSF;'
            GOTO 26
         ENDIF
*
*   Successfully read a CSF; update NREC and read another CSF
*
   17    NREC = NREC + 3
         GOTO 3
*
      ELSE  ! the record just read is either ' *' or EOF, marking
            ! the end of a block or end of the file
*
*   There is always at least one CSF
*
         IF (NCF .EQ. 1) THEN
            DO I = 1, NCORE
               CALL PACK (NKJ(I)+1, I, IQA(1,1))
               CALL PACK (0,        I, JQSA(1,1,1))
               CALL PACK (0,        I, JQSA(1,2,1))
               CALL PACK (1,        I, JQSA(1,3,1))
            ENDDO
            CALL PACK (0, 1,    JCUPA(1,1))
            CALL PACK (1, NNNW, JCUPA(1,1))
         ELSE
            NCF = NCF - 1
         ENDIF
*
      ENDIF

      IF (ncf .NE. ncfblock) THEN
         WRITE (istde,*) MYNAME//': ncf=',ncf, 'ncfblock=',ncfblock
         STOP
      ENDIF
*
*   Check if any subshell is empty; eliminate it from the
*   list if this is the case; issue a message
*
      I = NCORP1
   19 IF (I .LE. NW) THEN
         DO J = 1, NCF
            IF (IQ (I,J) .NE. 0) GOTO 23
         ENDDO
         CALL CONVRT (NP(I), str, LENTH)
         PRINT *, 'Subshell '//str(1:LENTH)//NH(I)//' is empty'
     &,                     ' in all CSFs'
   23    I = I + 1
         GOTO 19
      ENDIF
*
*   Store the number of electrons in the COMMON variable
*   This will act as a check now - it's been determined in lodcsh
*
      NCOREL = 0
      DO I = 1, NCORE
         NCOREL = NCOREL + NKJ(I) + 1
      ENDDO
!      NELEC = NCOREL+NPEEL
      IF (ncorel + npeel .NE. nelec) THEN
         WRITE (istde,*) MYNAME//': nelec not equal to that in lodcsh'
         STOP
      ENDIF

      if(myid==0) PRINT *, 'There are ',ncf,' 
     1   relativistic CSFs... load complete;'

      RETURN
*
   26 BACKSPACE (nfile)
   27 BACKSPACE (nfile)
   28 BACKSPACE (nfile)
      WRITE (istde,*) ' CSF sequence number: ', ncf
      DO 29 I = 1, 3
         READ  (nfile,'(A)',ERR = 29,END = 29) str
         WRITE (istde,*) str(1:LEN_TRIM (str))
   29 CONTINUE
   30 CLOSE (nfile)

      STOP
      END
