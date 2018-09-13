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
    1 PRINT *, 'Full name of the second CSF file'
      READ (*,'(A)') FILNAM
*
      IF ( LEN_TRIM (FILNAM) .EQ. 0) GOTO 1
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
*      DO 3 I = NCORER+1,NW
*         IFULLI = NKJ(I)+1
*         DO 2 J = 1,NCF
*            IF (IQ (I,J) .NE. IFULLI) GOTO 4
*    2    CONTINUE
*         CALL CONVRT (NP(I),RECORD,LENTH)
*         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
*         PRINT *, ' in all CSFs; including this'
*         PRINT *, ' subshell in the core;'
*         NCORE = NCORE+1
*    3 CONTINUE
*
    4 RETURN
      END
