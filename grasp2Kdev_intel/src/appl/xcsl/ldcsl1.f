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
    1 PRINT *, 'Full name of the first CSF file'
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
*      DO 3 I = NCORE+1,NW
*         IFULLI = NKJ(I)+1
*         DO 2 J = 1,NCF
*            IF (IQ (I,J) .NE. IFULLI) GOTO 4
*    2    CONTINUE
*         CALL CONVRT (NP(I),RECORD,LENTH)
*         PRINT *, 'Subshell '//RECORD(1:LENTH)//NH(I)//' is full'
*         PRINT *, ' in all CSFs; including this'
*         PRINT *, ' subshell in the core;'
*         NCORER = NCORER+1
*    3 CONTINUE
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
