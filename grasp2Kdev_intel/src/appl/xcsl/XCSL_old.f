************************************************************************
*                                                                      *
      SUBROUTINE CHKPLT
*                                                                      *
*   This code checks for  consistent substitutions of plants between   *
*   XCSL and the LIB92 subprograms.                                    *
*                                                                      *
*   Call(s) to: [LIB92]: CONVRT, LODPLT.                               *
*                                                                      *
*   Written by Farid A Parpia               Last update: 25 Sep 1993   *
*                                                                      *
************************************************************************
*
      LOGICAL LPLANT
      CHARACTER*256 RECORD
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
      COMMON/LIB92P/LPLANT,NPLANT(4)
*
*   Load COMMON/LIB92P/
*
      CALL LODPLT
*
*   Consistent numerical plants?
*
      IF (NPLANT(3) .NE. NNNW) THEN
         CALL CONVRT (NPLANT(3),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NW has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNW'
         PRINT *, ' in XCSL.'
         STOP
      ENDIF
*
      IF (NPLANT(4) .NE. NNNWP) THEN
         CALL CONVRT (NPLANT(4),RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NWP has been set to '
         PRINT *, ' '//RECORD(1:LENTH)//' in LIB92, but to NNNWP'
         PRINT *, ' in XCSL.'
         STOP
      ENDIF
*
      IF (MOD (NNNW,4) .EQ. 0) THEN
         NWP = NNNW/4
      ELSE
         NWP = NNNW/4+1
      ENDIF
*
      IF (NNNWP .NE. NWP) THEN
         CALL CONVRT (NWP,RECORD,LENTH)
         PRINT *, 'CHKPLT: Plant NWP should be set'
         PRINT *, ' to '//RECORD(1:LENTH)//'.'
         STOP
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      FUNCTION IQR (ISUBSH,ICSF)
*                                                                      *
*   IQR is the occupation of subshell ISUBSH in CSF  ICSF.             *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PNTIQR,IQAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
*
      IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NWR)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
            IQR = IUNPCK (IQAR(1,ICSF),ISUBSH)
         ELSE
            PRINT *, 'IQR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         PRINT *, 'IQR: Argument ISUBSH is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      FUNCTION ISPARR (ICSF)
*                                                                      *
*   ISPARR is the value of P for CSF number ICSF.                      *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
         ISPARR = IUNPCK (JCUPAR(1,ICSF),NNNW)
         ISPARR = SIGN (1,ISPARR)
      ELSE
         PRINT *, 'ISPARR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      FUNCTION ITJPOR (ICSF)
*                                                                      *
*   ITJPOR is the value of 2J+1 for CSF number ICSF.                   *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
         ITJPOR = IUNPCK (JCUPAR(1,ICSF),NNNW)
         ITJPOR = ABS (ITJPOR)
      ELSE
         PRINT *, 'ITJPOR: Argument ICSF is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      FUNCTION JCUPR (LOC,ICSF)
*                                                                      *
*   JCUPR is the 2J+1 value of the LOCth nontrivial intermediate ang-  *
*   ular momentum in CSF  ICSF.                                        *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNJQSR,PNTIQR
      POINTER (PNJQSR,JQSRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((LOC .GE. 1) .AND. (LOC .LE. NWR-1)) THEN
         IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
            JCUPR = IUNPCK (JCUPAR(1,ICSF),LOC)
         ELSE
            PRINT *, 'JCUPR: Argument ICSF is out of range.'
            STOP
         ENDIF
      ELSE
         PRINT *, 'JCUPR: Argument LOC is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      FUNCTION JQSR (IWHICH,ISUBSH,ICSF)
*                                                                      *
*   JQSR is a subshell quantum number for subshell ISUBSH in configu-  *
*   ration state function  ICSF:  the seniority if IWHICH is 1;  the   *
*   quantum number w if IWHICH is 2, and 2J+1 if IWHICH is 3.          *
*                                                                      *
*   Call(s) to: [LIB92]: IUNPCK.                                       *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PJCUPR,PNTIQR
      POINTER (PJCUPR,JCUPRDUMMY)
      POINTER (PNTIQR,IQRDUMMY)
*
      include 'parameters.def'
CGG      INTEGER NNNWP
CGG      PARAMETER (NNNWP = 30)
*
      POINTER (PNJQSR,JQSAR(NNNWP,3,*))
*
      COMMON/ORB2R/NCFR,NWR,PNTIQR
     :      /STATR/PNJQSR,PJCUPR
*
      IF ((IWHICH .GE. 1) .AND. (IWHICH .LE. 3)) THEN
         IF ((ISUBSH .GE. 1) .AND. (ISUBSH .LE. NWR)) THEN
            IF ((ICSF .GE. 1) .AND. (ICSF .LE. NCFR)) THEN
               JQSR = IUNPCK (JQSAR(1,IWHICH,ICSF),ISUBSH)
            ELSE
               PRINT *, 'JQSR: Argument ICSF is out of range.'
               STOP
            ENDIF
         ELSE
            PRINT *, 'JQSR: Argument ISUBSH is out of range.'
            STOP
         ENDIF
      ELSE
         PRINT *, 'JQSR: Argument IWHICH is out of range.'
         STOP
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE LDCSL1 (NCORER)
*                                                                      *
*   Open, check, load data from and close the  .csl  file. This file   *
*   is always attached to stream 21.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CGG      PARAMETER (NNNW = 120)
CGG      PARAMETER (NNNWP = 30)
      CHARACTER*256 FILNAM
      CHARACTER*15 RECORD
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
      CHARACTER*2 NH
*
      POINTER (PNTRIQ,IQA(NNNWP,*))
      POINTER (PNTIQR,IQAR(NNNWP,*))
      POINTER (PNTJQS,JQSA(NNNWP,3,*))
      POINTER (PNJCUP,JCUPA(NNNWP,*))
      POINTER (PNJQSR,JQSAR(NNNWP,3,*))
      POINTER (PJCUPR,JCUPAR(NNNWP,*))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF1R/NELECR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB2R/NCFR,NWR,PNTIQR
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB4R/NPR(NNNW),NAKR(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /STAT/PNTJQS,PNJCUP
     :      /STATR/PNJQSR,PJCUPR
*
*   The  .csl  file is FORMATTED; it must exist
*
      FORM = 'FORMATTED'
      STATUS = 'OLD'
*
*   Determine the name of the first  .csl  file
*
    1 PRINT *, 'Enter the name of the first GRASP92'
      PRINT *, ' Configuration Symmetry List File:'
      READ (*,'(A)') FILNAM
*
      IF (LENGTH (FILNAM) .EQ. 0) GOTO 1
*
      CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) GOTO 1
*
*   Check the first record of the file; if not as expected, try again
*
      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.
     :    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         PRINT *, 'Not a Configuration Symmetry List File;'
         CLOSE (21)
         GOTO 1
      ENDIF
*
*   Load data from the  .csl  file
*
      CALL LODCSL (NCORE)
*
*   Close the  .csl  file
*
      CLOSE (21)
*
*   Check if the core should be redefined
*
      NCORER = NCORE
      DO 3 I = NCORE+1,NW
         IFULLI = NKJ(I)+1
         DO 2 J = 1,NCF
            IF (IQ (I,J) .NE. IFULLI) GOTO 4
    2    CONTINUE
         CALL CONVRT (NP(I),RECORD,LENTH)
         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
         PRINT *, ' in all CSFs; including this'
         PRINT *, ' subshell in the core;'
         NCORER = NCORER+1
    3 CONTINUE
*
*   Copy data to reference arrays
*
    4 NELECR = NELEC
      NCFR = NCF
      NWR = NW
*
      CALL ALLOC (PNTIQR,NNNWP  *NCF,4)
      CALL ALLOC (PNJQSR,NNNWP*3*NCF,4)
      CALL ALLOC (PJCUPR,NNNWP  *NCF,4)
*
      DO 6 I = 1,NCF
         DO 5 J = 1,NNNWP
            IQAR(J,I) = IQA(J,I)
            JQSAR(J,1,I) = JQSA(J,1,I)
            JQSAR(J,2,I) = JQSA(J,2,I)
            JQSAR(J,3,I) = JQSA(J,3,I)
            JCUPAR(J,I) = JCUPA(J,I)
    5    CONTINUE
    6 CONTINUE
*
      DO 7 I = 1,NW
         NPR(I) = NP(I)
         NAKR(I) = NAK(I)
    7 CONTINUE
*
*   Deallocate storage so that LODCSL can restart correctly
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE LDCSL2 (NCORE)
*                                                                      *
*   Open, check, load data from and close the  .csl  file. This file   *
*   is always attached to stream 21.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: IQ, LENGTH, LODCSL, OPENFL.                   *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
Cww      INTEGER PNTRIQ
      POINTER (PNTRIQ,RIQDUMMY)
      CHARACTER*256 FILNAM
      CHARACTER*15 RECORD
      CHARACTER*11 FORM
      CHARACTER*3 STATUS
      CHARACTER*2 NH
*
      include 'parameters.def'
CGG      INTEGER NNNW
CGG      PARAMETER (NNNW = 120)
      COMMON/ORB2/NCF,NW,PNTRIQ
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
*
*   The  .csl  file is FORMATTED; it must exist
*
      FORM = 'FORMATTED'
      STATUS = 'OLD'
*
*   Determine the name of the first  .csl  file
*
    1 PRINT *, 'Enter the name of the second GRASP92'
      PRINT *, ' Configuration Symmetry List File:'
      READ (*,'(A)') FILNAM
*
      IF (LENGTH (FILNAM) .EQ. 0) GOTO 1
*
      CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) GOTO 1
*
*   Check the first record of the file; if not as expected, try again
*
      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.
     :    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         PRINT *, 'Not a Configuration Symmetry List File;'
         CLOSE (21)
         GOTO 1
      ENDIF
*
*   Load data from the  .csl  file
*
      CALL LODCSL (NCORER)
*
*   Close the  .csl  file
*
      CLOSE (21)
*
*   Check if the core should be redefined
*
      NCORE = NCORER
      DO 3 I = NCORER+1,NW
         IFULLI = NKJ(I)+1
         DO 2 J = 1,NCF
            IF (IQ (I,J) .NE. IFULLI) GOTO 4
    2    CONTINUE
         CALL CONVRT (NP(I),RECORD,LENTH)
         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
         PRINT *, ' in all CSFs; including this'
         PRINT *, ' subshell in the core;'
         NCORE = NCORE+1
    3 CONTINUE
*
    4 RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE SETDBG
*                                                                      *
*   This subroutine sets the arrays that control debug printout from   *
*   the radial and angular modules of the GRASP92 suite.               *
*                                                                      *
*   Call(s) to: [LIB92]: GETYN, LENGTH, OPENFL.                        *
*                                                                      *
*   Written by Farid A Parpia               Last update: 23 Dec 1992   *
*                                                                      *
************************************************************************
*
      LOGICAL GETYN,LDBPA,LDBPG,LDBPR,YES
      CHARACTER*256 FILNAM
      CHARACTER*11 DEFNAM,FORM
      CHARACTER*3 STATUS
*
      COMMON/DEBUGA/LDBPA(5)
     :      /DEBUGG/LDBPG(5)
     :      /DEBUGR/LDBPR(30)
*
*   Initialise the arrays that control the debug printout
*
      DO 1 I = 1,5
         LDBPA(I) = .FALSE.
    1 CONTINUE
*
      DO 2 I = 1,5
         LDBPG(I) = .FALSE.
    2 CONTINUE
*
      DO 3 I = 1,30
         LDBPR(I) = .FALSE.
    3 CONTINUE
*
      PRINT *, 'Generate debug printout?'
      YES = GETYN ()
      IF (YES) THEN
*
*   The  .dbg  file is formatted; open it on unit 99
*
         DEFNAM = 'rci92.dbg'
         FORM = 'FORMATTED'
         STATUS = 'NEW'
*
         PRINT *, 'File  rci92.dbg  will be created as the'
         PRINT *, ' RCI92 DeBuG Printout File; enter another'
         PRINT *, ' file name if this is not acceptable;'
         PRINT *, ' null otherwise:'
         READ (*,'(A)') FILNAM
*
         IF (LENGTH (FILNAM) .EQ. 0) FILNAM = DEFNAM
*
    4    CALL OPENFL (99,FILNAM,FORM,STATUS,IERR)
         IF (IERR .NE. 0) THEN
    5       PRINT *, 'Enter a name for the RCI92 DeBuG Printout'
            PRINT *, ' file that is to be created:'
            READ (*,'(A)') FILNAM
            IF (LENGTH (FILNAM) .EQ. 0) GOTO 5
            GOTO 4
         ENDIF
*
*   Set options for general printout
*
         PRINT *, ' Print out the machine constants used?'
         YES = GETYN ()
         IF (YES) LDBPG(1) = .TRUE.
*
*   Set options for angular modules
*
         PRINT *, ' Printout from LODCSL?'
         YES = GETYN ()
         IF (YES) LDBPA(1) = .TRUE.
*
      ENDIF
*
      RETURN
      END
************************************************************************
*                                                                      *
      SUBROUTINE XCLD12 (NCORER,NCORE)
*                                                                      *
*   Compares the second .csl list with the first .csl list. In addi-   *
*   tion to all  CSFs in the reference list,  all  CSFs that are not   *
*   repeated in the second .csl list are written out.                  *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, CONVRT, IQ, ISPAR, ITJPO, JCUP, JQS,   *
*                        LENGTH, OPENFL.                               *
*               [XCSL]: IQR, ISPARR, ITJPOR, JCUPR, JQSR, OUTCSF.      *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 28 Sep 1993   *
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
      CHARACTER*11 CNUMBR,FORM
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
*   The reference  .csl  list must have more or an equal number of
*   subshells
*
      IF (NWR .LT. NW) THEN
         PRINT *, 'Too many subshells in the second GRASP92'
         PRINT *, ' Configuration Symmetry List File: this'
         PRINT *, ' should be the first File.'
         STOP
      ENDIF
*
*   Check to see that the subshell symmetry lists are ordered
*   identically as far as possible
*
      DO 1 I = 1,NW
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
    2 PRINT *, 'Enter the name of the GRASP92'
      PRINT *, ' Configuration Symmetry List'
      PRINT *, ' File that is to be created:'
      READ (*,'(A)') FILNAM
*
      IF (LENGTH (FILNAM) .EQ. 0) GOTO 2
*
      CALL OPENFL (21,FILNAM,FORM,STATUS,IERR)
      IF (IERR .EQ. 1) GOTO 2
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
      IF ((NWR - NCORE).LE.80) THEN
         WRITE (21,300) (NP(I),NH(I),I = NCORE+1,NWR)
      ELSE
         WRITE (21,300) (NP(I),NH(I),I = NCORE+1,NCORE+80)
         WRITE (21,300) (NP(I),NH(I),I = NCORE+81,NWR)
      ENDIF
      WRITE (21,'(A)') 'CSF(s):'
*
*   Write out all CSFs which occur in the reference list but not in
*   the second list
*
      ICOUNT = 0
*
      DO 5 I = 1,NCFR
        ITJPI = ITJPOR (I)
        ISPAI = ISPARR (I)
        DO 4 J = 1,NCF
          IF (ITJPI .NE. ITJPO (J)) GOTO 4
          IF (ISPAI .NE. ISPAR (J)) GOTO 4
          NOPEN = 0
          DO 3 K = 1,NW
            IQKJ = IQ (K,J)
            IF (IQR (K,I) .NE. IQKJ) GOTO 4
            IF (JQSR (1,K,I) .NE. JQS (1,K,J)) GOTO 4
            IF (JQSR (2,K,I) .NE. JQS (2,K,J)) GOTO 4
            IF (JQSR (3,K,I) .NE. JQS (3,K,J)) GOTO 4
            IF ((IQKJ .GT. 0) .AND. (IQKJ .LT. NKJ(K)+1))
     :        NOPEN = NOPEN+1
            IF (NOPEN .GT. 1) THEN
               LOC = NOPEN-1
               IF (JCUPR (LOC,I) .NE. JCUP (LOC,J)) GOTO 4
            ENDIF
    3     CONTINUE
          GOTO 5
    4   CONTINUE
        ICOUNT = ICOUNT+1
        CALL OUTCSF (I,NCORE,NWR,IQR,ISPARR,ITJPOR,JCUPR,JQSR)
    5 CONTINUE
*
      CALL CONVRT (NWR,CNUMBR,LENTH)
      PRINT *, CNUMBR(1:LENTH)//' relativistic subshells;'
      CALL CONVRT (ICOUNT,CNUMBR,LENTH)
      PRINT *, CNUMBR(1:LENTH)//' relativistic CSFs.'
*
      CLOSE (21)
*
      RETURN
*
  300 FORMAT (51(1X,1I2,1A2))
*
      END
************************************************************************
************************************************************************
************************************************************************
***                                                                  ***
***                                                                  ***
***                **   **   *****    *****   **                     ***
***                 ** **   **   **  **   **  **                     ***
***                  ***    **       **       **                     ***
***                  ***    **        *****   **                     ***
***                  ***    **            **  **                     ***
***                 ** **   **   **  **   **  **                     ***
***                **   **   *****    *****   *******                ***
***                                                                  ***
***   Utility to delete a set of CSFs from the CSL File              ***
***                                                                  ***
***                            GRASP92                               ***
***          F. A. Parpia, C. F. Fischer, and I. P. Grant            ***
***                                                                  ***
************************************************************************
************************************************************************
************************************************************************
*                                                                      *
      PROGRAM XCSL
*                                                                      *
*   Entry routine for XCSL. Controls the entire computation.           *
*                                                                      *
*   Call(s) to: [LIB92]: SETMC.                                        *
*               [XCSL]: CHKPLT, LDCSL1, LDCSL2, XCLD12, SETDBG.        *
*                                                                      *
*   Written by Farid A. Parpia            Last revision: 25 Sep 1993   *
*                                                                      *
************************************************************************
*
      PRINT *, 'XCSL: Execution begins ...'
*
*   Check compatibility of plant substitutions
*
      CALL CHKPLT
*
*   Determine if there is to be any debug printout; this will be
*   made on the  .dbg  file
*
      CALL SETDBG
*
*   Perform machine- and installation-dependent setup
*
      CALL SETMC
*
*   Load the first  .csl  file
*
      CALL LDCSL1 (NCORER)
*
*   Load the second  .csl  file
*
      CALL LDCSL2 (NCORE)
*
*   Eliminate CSFs common to both lists
*
      CALL XCLD12 (NCORER,NCORE)
*
*   Print completion message
*
      PRINT *, 'XCSL: Execution complete.'
*
      STOP
      END
