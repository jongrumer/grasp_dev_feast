************************************************************************
*                                                                      *
      SUBROUTINE MERG12 (NCORER,NCORE)
*                                                                      *
*   Compares the second .csl list with the first .csl list. In addi-   *
*   tion to all  CSFs in the reference list,  all  CSFs that are not   *
*   repeated in the second .csl list are written out.                  *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
*                        LENGTH, OPENFL.                               *
*               [MRGCSL]: IQR, ISPARR, ITJPOR, JCUPR, JQSR, OUTCSF.    *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 24 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
Cww      INTEGER PNTIQR,PNTRIQ
      POINTER (PNTIQR,IQRDUMMY)
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 FILNAM
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
      CHARACTER*2 NH
*
      EXTERNAL IQ,IQR,ISPAR,ISPARR,ITJPO,ITJPOR,JCUP,JCUPR,JQS,JQSR
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF1R/NELECR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB2R/NCFR,NWR,PNTIQR
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB4R/NPR(NNNW),NAKR(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
*   The same number of electrons must appear in both lists
*
      IF (NELECR .NE. NELEC) THEN
         PRINT *, 'The number of electrons is not equal in the'
         PRINT *, ' first and second GRASP92 Configuration'
         PRINT *, ' symmetry list files.'
         STOP
      ENDIF
*
*   The reference  .csl  list must have fewer or an equal number of
*   subshells
*
      IF (NWR .GT. NW) THEN
         PRINT *, 'Too many subshells in the second GRASP92'
         PRINT *, ' Configuration Symmetry List File: this'
         PRINT *, ' should be the first File.'
         STOP
      ENDIF
*
*   Check to see that the subshell symmetry lists are ordered
*   identically as far as possible
*
      DO 1 I = 1,NWR
         IF ((NPR(I) .NE. NP(I)) .OR. (NAKR(I) .NE. NAK(I))) THEN
            PRINT *, 'Subshells must appear in the same order in both'
            PRINT *, ' GRASP92 Configuration Symmetry List Files.'
            STOP
         ENDIF
    1 CONTINUE
*
*   Open the output  .csl file; it is FORMATTED; it must not exist
*
      FORM = 'FORMATTED'
      STATUS = 'NEW'
*
*   Determine the name of the output  .csl  file
*
C    2 PRINT *, 'Name of the CSF file to be created'
      FILNAM = 'rcsf.out'
*
C      IF ( LEN_TRIM (FILNAM) .EQ. 0) GOTO 2
*
      CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
C      IF (IERR .EQ. 1) GOTO 2
*
*   Determine the common core subshells; write out the list;
*   determine the pell subshells; write out the list; these
*   data form the first part of the header of the .csl file;
*   one additional line forms the remainder of the header of
*   the  .csl file
*
      NCORE = MIN (NCORE,NCORER)
      WRITE (21,'(A)') 'Core subshells:'
      WRITE (21,300) (NP(I),NH(I),I = 1,NCORE)
      WRITE (21,'(A)') 'Peel subshells:'
      IF ((NW - NCORE).LE.80) THEN
         WRITE (21,300) (NP(I),NH(I),I = NCORE+1,NW)
      ELSE
         WRITE (21,300) (NP(I),NH(I),I = NCORE+1,NCORE+80)
         WRITE (21,300) (NP(I),NH(I),I = NCORE+81,NW)
      ENDIF
      WRITE (21,'(A)') 'CSF(s):'
*
*   Now write out all CSFs in the reference list
*
      DO 4 I = 1,NCFR
         CALL OUTCSF (I,NCORE,NWR,IQR,ISPARR,ITJPOR,JCUPR,JQSR)
    4 CONTINUE
*
*   Now write out all unrepeated CSFs in the second list
*
      ICOUNT = NCFR
*
      DO 7 I = 1,NCF
*
         DO 6 J = 1,NCFR
            NOPEN = 0
            DO 5 K = 1,NWR
               IQRKJ = IQR (K,J)
               IF ((IQRKJ .GT. 0) .AND. (IQRKJ .LT. NKJ(K)+1))
     :            NOPEN = NOPEN+1
               IF ((IQRKJ .NE. IQ (K,I)) .OR.
     :             (JQSR (1,K,J) .NE. JQS (1,K,I)) .OR.
     :             (JQSR (2,K,J) .NE. JQS (2,K,I)) .OR.
     :             (JQSR (3,K,J) .NE. JQS (3,K,I))) THEN
                  GOTO 6
               ELSE
                  IF (NOPEN .GT. 1) THEN
                     LOC = NOPEN-1
                     IF (JCUPR (LOC,J) .NE. JCUP (LOC,I)) GOTO 6
                  ENDIF
               ENDIF
    5       CONTINUE
            GOTO 7
    6    CONTINUE
*
         ICOUNT = ICOUNT+1
         CALL OUTCSF (I,NCORE,NW,IQ,ISPAR,ITJPO,JCUP,JQS)
    7 CONTINUE
*
      CALL CONVRT (NW,FILNAM,LENTH)
      PRINT *, FILNAM(1:LENTH)//' relativistic subshells;'
      CALL CONVRT (ICOUNT,FILNAM,LENTH)
      PRINT *, FILNAM(1:LENTH)//' relativistic CSFs.'
*
      CLOSE (21)
*
      RETURN
*
  300 FORMAT (120(1X,1I2,1A2))
*
      END
