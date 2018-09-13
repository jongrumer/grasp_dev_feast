************************************************************************
*                                                                      *
      SUBROUTINE LDCSL1 (NCORER,NAME)
*                                                                      *
*   Open, check, load data from and close the  .csl  file. This file   *
*   is always attached to stream 21.                                   *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, IQ, LENGTH, LODCSL, OPENFL.            *
*                                                                      *
*                                                                      *
************************************************************************
*
      IMPLICIT DOUBLEPRECISION (A-H, O-Z)
      include 'parameters.def'
CFF      PARAMETER (NNNW = 120)
CFF      PARAMETER (NNNWP = 30)
      CHARACTER*24 NAME
      CHARACTER*15 RECORD
      CHARACTER*2 NH, NHR, NHII
*
      INTEGER*4 IQA, JQSA, JCUPA
      INTEGER*4 IQAR,JQSAR,JCUPAR
 
      POINTER (PNTRIQ,IQA(NNNWP,1))
      POINTER (PNTIQR,IQAR(NNNWP,1))
      POINTER (PNTJQS,JQSA(NNNWP,3,1))
      POINTER (PNJCUP,JCUPA(NNNWP,1))
      POINTER (PNJQSR,JQSAR(NNNWP,3,1))
      POINTER (PJCUPR,JCUPAR(NNNWP,1))
*
      COMMON/DEF1/EMN,IONCTY,NELEC,Z
     :      /DEF1R/NELECR
     :      /ORB2/NCF,NW,PNTRIQ
     :      /ORB2R/NCFR,NWR,PNTIQR
     :      /ORB4/NP(NNNW),NAK(NNNW)
     :      /ORB4R/NPR(NNNW),NAKR(NNNW)
     :      /ORB5/NKL(NNNW),NKJ(NNNW)
     :      /ORB10/NH(NNNW)
     :      /ORB10R/NHR(NNNW)
     :      /STAT/PNTJQS,PNJCUP
     :      /STATR/PNJQSR,PJCUPR

*
*     Common relevant for the radial wave functions of the initial state.
*
      COMMON/DEF1II/EMNII,IONCTYII,NELECII,ZII
      COMMON/ORB2II/NCFII,NWII
     :      /ORB4II/NPII(NNNW),NAKII(NNNW)
     :      /ORB5II/NKLII(NNNW),NKJII(NNNW)
     :      /ORB10II/NHII(NNNW)
*
*   Check the first record of the file; if not as expected, try again
*
      J = INDEX(NAME,' ')
      OPEN (UNIT = 21,FILE=NAME(1:J-1)//'.c',FORM='FORMATTED',
     :     STATUS='OLD')

      READ (21,'(1A15)',IOSTAT = IOS) RECORD
      IF ((IOS .NE. 0) .OR.
     :    (RECORD(1:15) .NE. 'Core subshells:')) THEN
         PRINT *, 'Not a Configuration Symmetry List File;'
         CLOSE (21)
         STOP
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
      NELECII = NELEC
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
      NWII=NW
      NCFII=NCF
      DO 7 I = 1,NW
         NHII(I) = NH(I)
         NPII(I) = NP(I)
         NAKII(I) = NAK(I)
         NKLII(I) = NKL(I)
         NKJII(I) = NKJ(I)
         NPR(I) = NP(I)
         NAKR(I) = NAK(I)
         NHR(I) = NH(I)
    7 CONTINUE
*
*   Deallocate storage so that LODCSL can restart correctly
*
      CALL DALLOC (PNTRIQ)
      CALL DALLOC (PNTJQS)
      CALL DALLOC (PNJCUP)
*
      CLOSE (21)

      RETURN
      END
